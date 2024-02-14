# Author: Max Freudenberg
# maximilian.freudenberg@uni-goettingen.de

# The purpose of this script is to extract sentinel pixel time series
# from a set of pre-cutout images. The pixels are cut out at the NFI
# tree positions. Each tree comes with a modeled disk, approximating 
# the size it covers. The final pixel values per disk are calculated 
# as the weighted mean of the intersections with the underlying pixels.

# How it works:
# 1) Load all trees that existed in 2012 and still exist in 2022 from
#    external postgres db.
# 2) Add non-tree positions to be cut out as well (by closest point)
# 3) Load all images belonging to a certain NFI plot
# 4) Load all crown disk polygons belonging to the plot
# 5) Iterate over trees / polygons and their underlying pixels,
#    cut them out and calculate weighted average
# 6) Add information on stand purity
# 7) Add train / test split
# 8) Save

using Base.Iterators
using ArgParse
using Statistics
using Serialization
using ForceCubeAccess
using Dates
using CSV
using Base.Threads
using Base.Threads: atomic_add!
import GeoInterface as GI
import ArchGDAL as AG
using Rasters
using DimensionalData
include("utils.jl")
include("pixel_operations.jl")

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--output", "-o"
            help = "Final database output file path including .sqlite extension."
            required = true
        "--tilling", "-t"
            help = "Path to CSV file containing tilling data (Bestockung) for each plot."
            required = true
        "--nonforest-points", "-n"
            help = "Path to sqlite file containing geolocations of non-forest points."
            required = true
        "--tile-folder", "-f"
            help = "Path to folder containing Sentinel-2 cutout subfolders. Each cutout 
            must be a subfolder named by the cluster plot id, containing the tiles."
            required = true
        "--pg-database", "-d"
            help = "Connection to postgres database containing NFI tree information. 
            Has to be specified in the form postgresql://username:pw@ip:port/databasename"
            required = true
    end

    return parse_args(s)
end

"""
    extract_pixels_worker(path, chunk, i; lazy = false)

The main pixel extraction worker routine. Acts on chunks of an input DataFrame 
containing the tree id, species and coordinates.

Args:
- path: Path to the folder containing the S2 BOA and QAI tiles
- chunk: DataFrame grouped by cluster plot id (tnr) with columns tree_id, ba, 
    geom and geom_sfl.
- i: Atomic integer to keep track of progress in threaded computation.
- lazy: Whether to open the underlying rasters lazily.
"""
function extract_pixels_worker(path, chunk, i; lazy = false)
    res = DataFrame(:tnr => Int32[],
        :enr => Int32[],
        :tree_id => Int32[],
        :ba => Int32[],
        :time => DateTime[],
        :boa => Vector{Int16}[],
        :qai => UInt16[])
    missed_plots = 0
    missed_trees = 0
    for group in chunk
        tnr = first(group.tnr)
        # @info tnr
        local boas
        local qais
        local boa_series
        local qai_series
        try
            boas = filter(contains("BOA"),
                readdir(joinpath(path, string(tnr)); join = true))
            qais = filter(contains("QAI"),
                readdir(joinpath(path, string(tnr)); join = true))
            times = DateTime.(getindex.(split.(basename.(boas), '_'), 2))
            boa_series = RasterSeries(boas, Ti(times); duplicate_first = true, lazy = lazy)
            qai_series = RasterSeries(qais, Ti(times); duplicate_first = true, lazy = lazy)
        catch err
            if err isa Base.IOError
                missed_plots += 1
                missed_trees += nrow(group)
                atomic_add!(i, 1)
                continue
            elseif err isa DimensionalData.DimensionMismatch
                @warn "Dimension mismatch on plot $tnr"
                atomic_add!(i, 1)
                continue
            else
                rethrow(err)
            end
        end

        sr = boa_series[1]  # sample raster
        xmin_global, ymin_global = minimum.(dims(sr, (X, Y)))
        xstep, ystep = step.(dims(sr, (X, Y)))

        corner_groups = groupby(group, :enr)

        for corner in corner_groups
            enr = first(corner.enr)
            # non-tree points have a corner id (enr) of -1
            if enr == -1
                for (point, tree_id, ba) in zip(corner.geom,
                    corner.tree_id,
                    corner.ba)
                    x, y = GI.coordinates(point)
                    extracted_boas = map(boa_series) do raster
                        return Float32.(raster[X(Near(x)), Y(Near(y))])
                    end
                    extracted_qais = map(qai_series) do raster
                        return raster[X(Near(x)), Y(Near(y))]
                    end
                    for (time, boa, qai) in zip(dims(extracted_boas, Ti),
                        extracted_boas,
                        extracted_qais)
                        push!(res, [tnr, enr, tree_id, ba, time, boa, qai])
                    end
                end
                # tree extraction
            else
                valid_polygons = valid_polygon_indices(corner.geom_sfl, 0.5, 3)
                corner = corner[valid_polygons, :]

                for (polygon, tree_id, ba) in zip(corner.geom_sfl,
                    corner.tree_id,
                    corner.ba)
                    grid = underlying_pixelgrid(polygon, xmin_global, ymin_global, xstep,
                        ystep)
                    pixel_polygons = coord_to_polygon.(grid, xstep, ystep)
                    extracted_boas = map(boa_series) do raster
                        return weighted_mean_pixel_value(polygon,
                            grid,
                            pixel_polygons,
                            raster)
                    end
                    # just take the closest qai value as proxy
                    x, y = grid[1]
                    extracted_qais = map(qai_series) do raster
                        return raster[X(Near(x)), Y(Near(y))]
                    end
                    # return extracted_timeseries
                    for (time, boa, qai) in zip(dims(extracted_boas, Ti),
                        extracted_boas,
                        extracted_qais)
                        push!(res, [tnr, enr, tree_id, ba, time, boa, qai])
                    end
                end
            end
        end
        atomic_add!(i, 1)
        if i[] % 500 == 0
            GC.gc()
        end
    end
    return res, missed_plots, missed_trees
end

"""
    extract_pixels(path, df)

Parallelized pixel extraction routine. Needs the path to the S2 BOA and QAI tiles 
as well as a DataFrame with columns tree_id, ba, geom, geom_sfl.
"""
function extract_pixels(path, df)
    chunksize = length(df) < Threads.nthreads() ? 1 : length(df) ÷ Threads.nthreads() ÷ 4
    chunks = Iterators.partition(df, chunksize)
    i = Threads.Atomic{Int}(0)
    last_i = 0

    tasks = map(chunks) do chunk
        Threads.@spawn extract_pixels_worker(path, chunk, i)
    end

    while i[] < length(df)
        sleep(5)
        increment = i[] - last_i
        last_i = i[]
        println("$(i[]) / $(length(df)) \t $(increment/5)it/s")
    end

    results = fetch.(tasks)
    missed_plots = sum(getindex.(results, 2))
    missed_trees = sum(getindex.(results, 3))
    @warn "Missed $missed_plots plots and $missed_trees trees."
    return vcat(first.(results)...)
end

"""
    add_purity_info!(df, csv_file)

Reads a CSV containing info on whether a subplot is pure or not and 
attaches this info on DataFrame df. `tree_data` is the raw DataFrame 
containing the tree positions and subplot ids.
"""
function add_purity_info!(df, csv_file)
    bestockung = CSV.read(csv_file, DataFrame)
    rename!(bestockung, :Tnr => :tnr)
    rename!(bestockung, :Enr => :enr)
    replace!(bestockung.BestockTypFein, "NULL" => "-1")
    bestockung.type = parse.(Int, bestockung.BestockTypFein)
    select!(bestockung, [:tnr, :enr, :type])
    bestockung_grouped = groupby(bestockung, [:tnr, :enr])

    df.is_pure = zeros(Bool, nrow(df))
    df.is_pure[df.ba .== -1, :] .= true

    df_grouped = groupby(df, [:tnr, :enr])

    for grp in df_grouped
        tnr, enr = first(grp.tnr), first(grp.enr)
        # skip if enr is negative (non-tree)
        enr < 0 && continue
        is_pure = first(bestockung_grouped[(tnr = tnr, enr = enr)].type) % 100 == 0
        grp.is_pure[:] .= is_pure
    end
end

"""
    train_test_partition!(df::AbstractDataFrame, train_split = 0.7)

Adds a column is_train to the DataFrame and sets values to one or zero 
according to the training split. The DataFrame has to have the cluster 
plot id as column with name tnr - records are assigned to train or test 
based on cluster plot id.
"""
function train_test_partition!(df::AbstractDataFrame; train_split = 0.7)
    df.is_train = fill(Bool(0), nrow(df))
    tnr_groups = groupby(df, :tnr)
    is_train = rand(length(tnr_groups)) .< train_split

    for (i, g) in enumerate(tnr_groups)
        g.is_train[:] .= is_train[i]
    end
end

"""
    main()

Generates the training dataset.
"""
function main()
    # args = parse_commandline()

    outfile = "dataset_V0.sqlite"
    nonforest_points = "datasets/nonforest_pixels.sqlite"
    s2_cutout_path = "/data_local_ssd/bwi/cutouts_2023/"
    purity_info_csv = "/data_hdd/bwi/bwi2012_bestockung.csv"

    # outfile = args["output"]
    # nonforest_points = args["nonforest-points"]
    # s2_cutout_path = args["tile-folder"]
    # purity_info_csv = args["tilling"]
    # db = args["pg-database"]
    db = ENV["PG_DATABASE"]

    if isfile(outfile)
        error("Output file $outfile exists. Please provide a different output file name.")
    end

    println("Loading data from postgres")
    # ENV["PG_DATABASE"] = postgresql://username:pw@ip:port/databasename
    conn = AG.read(db;
        flags = AG.OF_READONLY | AG.OF_VERBOSE_ERROR | AG.OF_VECTOR)

    # takes ~100s
    tree_data_raw = executesql(conn,
        "select * from geo.persistent_trees_matview where bhd>200 and bs!=2 and hoehe>120")
    # throw out trees without attached geometry
    tree_data = tree_data_raw[.!ismissing.(tree_data_raw.geom), :]

    # widen the type of geom_sfl to include missings, that we add further down
    tree_data.geom_sfl = Vector{Union{eltype(tree_data.geom_sfl), Missing}}(tree_data.geom_sfl)
    unique!(tree_data, :tree_id)
    sort!(tree_data, :bnr)

    println("Finished loading tree table")

    # now add coordinates that are not trees
    non_tree_coords = DataFrame(AG.getlayer(AG.read(nonforest_points), 0))
    rename!(non_tree_coords, :plot_id => :tnr, :GEOMETRY => :geom)
    for r in eachrow(non_tree_coords)
        push!(tree_data,
            (r.geom, missing, r.tree_id, -1, r.tnr, -1, -1, -1, -1, -1, -1, "", "", 0))
    end
    println("Nonforest info added")

    # tree_data = @view tree_data[1:1000:end, :]

    # group df by cluster plot id (tnr) and extract the pixels
    grouped_df = groupby(tree_data, :tnr)
    @time res = extract_pixels(s2_cutout_path, grouped_df)

    println("Adding purity info")
    @time add_purity_info!(res, purity_info_csv)
    println("Splitting train / test")
    @time train_test_partition!(res; train_split = 0.7)

    select!(res, Not([:tnr, :enr]))

    prepare_df_dtypes_for_export!(res)
    println("Writing result")
    @time write_sqlite(outfile, res)
    println("Done")
end

main()
