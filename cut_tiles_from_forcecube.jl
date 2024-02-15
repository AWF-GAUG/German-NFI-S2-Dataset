# for each tile
# look which rasters fall within the tile
# look into which blocks they fall
# group all rasters within one block
# make a bounding box of all grouped rasters
# read the qai for the entire bounding box
# for each raster, check if there is valid data in the qai
# make a new bounding box with all valid rasters
# read the entire boa for the valid bounding box
# crop each raster
# write raster to respective folder
using Glob: GlobMatch
using ForceCubeAccess
import ForceCubeAccess as FCA
using ArchGDAL
import ArchGDAL as AG
using Rasters
using Extents
using GDAL
using Serialization
import DimensionalData as DD
using DataFrames
using Dates
using ChunkSplitters


const cloudy = CLOUD_BUFFER | CLOUD_CIRRUS | CLOUD_OPAQUE | CLOUD_SHADOW
const cloudy_or_nodata = cloudy | NODATA
const mymask = cloudy_or_nodata | ILLUMIN_LOW | ILLUMIN_NONE | ILLUMIN_POOR | SNOW


function reprojected_bounds(r::Raster, target_crs)
    ((xmin, xmax),(ymin, ymax)) = Rasters.bounds(r, (X,Y))
    pts = [(xmin, ymin),(xmax, ymin),(xmax, ymax),(xmin, ymax)]
    # reproject from raster crs to target crs
    bds = AG.reproject(pts, crs(r), target_crs; order=:trad)
    xs = first.(bds)
    ys = last.(bds)
    xmin_, xmax_ = extrema(xs)
    ymin_, ymax_ = extrema(ys)
    return ((xmin_, xmax_), (ymin_, ymax_))
end


function tile_ranges(r::Raster, fc::ForceCube; margin=30)
    ((xmin, xmax),(ymin, ymax)) = reprojected_bounds(r, crs(fc))
    xmin -= margin
    xmax += margin
    ymin -= margin
    ymax += margin
    pts = ((xmin, ymin),(xmax, ymin),(xmax, ymax),(xmin, ymax))
    indices = (ForceCubeAccess.tile_index(fc, point...) for point in pts)
    xs = extrema(getindex.(indices, 1))
    ys = extrema(getindex.(indices, 2))
    return range(xs...), range(ys...)
end


function block_ranges(small_raster, big_raster, blocksize=(3000, 300))
    ((xmin, xmax),(ymin, ymax)) = reprojected_bounds(small_raster, crs(big_raster))
    
    xpxlims = extrema(DD.Dimensions.dims2indices(dims(big_raster,X), X(xmin..xmax)))
    ypxlims = extrema(DD.Dimensions.dims2indices(dims(big_raster,Y), Y(ymin..ymax)))

    xblocklims = ceil.(Int, xpxlims ./ blocksize[1])
    yblocklims = ceil.(Int, ypxlims ./ blocksize[2])
    return xblocklims, yblocklims
end


function getraster(raster, fc)
    xr, yr = tile_ranges(raster, fc)
    return first(fc[yr[1], xr[1]])
end


function get_access_groups(fc, rasters)
    # calculate the tile(s) on which the raster falls
    # raster bounds are reprojected to match the fc crs
    # some margin is added to be safe after reprojection
    tr = [tile_ranges(r, fc; margin=30) for r in rasters]

    df = DataFrame(:fn => Rasters.filename.(rasters), :r => rasters, :xr => getindex.(tr, 1), :yr => getindex.(tr, 2));
    df.onone = length.(df.xr) .== 1 .& length.(df.yr) .== 1;
    df.tnr = first.(splitext.(basename.(df.fn)))
    df.extent = [reprojected_extent_with_margin(r, crs(fc)) for r in rasters]

    # all rasters that are on one cube tile
    onone = filter(df) do row
        row.onone
    end;
    
    # all rasters that cover more than one tile
    others = filter(df) do row
        !row.onone
    end

    # extract the force cube tile for each raster so that we can check the block accesses
    big_rasters = [getraster(r, fc) for r in onone.r];
    blockrange = block_ranges.(onone.r, big_rasters);

    onone.blockrange = blockrange#getindex.(br,2)

    # access_groups = groupby(onone, [:xr, :yr]);
    access_groups = onone

    return access_groups, others
end


function combined_extent(rasters, crs; margin=30)
    (xmin, xmax), (ymin, ymax) = foldl(Extents.union, (bounds2extent(reprojected_bounds(r, crs)) for r in rasters))
    xmin -= margin
    xmax += margin
    ymin -= margin
    ymax += margin
    return Extent(X=(xmin,xmax), Y=(ymin,ymax)) 
end


function extent_with_margin(raster, margin=10)
    (xmin, xmax), (ymin, ymax) = Extents.extent(raster)
    xmin -= margin
    xmax += margin
    ymin -= margin
    ymax += margin
    return Extent(X=(xmin,xmax), Y=(ymin,ymax))
end


function reprojected_extent_with_margin(raster, crs; margin=30)
    (xmin, xmax), (ymin, ymax) = bounds2extent(reprojected_bounds(raster, crs))
    xmin -= margin
    xmax += margin
    ymin -= margin
    ymax += margin
    return Extent(X=(xmin,xmax), Y=(ymin,ymax))
end


function eachband(r::Raster)
    bands = dims(r, Band)
    return (view(r, Band(At(b))) for b in bands)
end


xyextent(ex::Extent) = Extent((X=ex.X, Y=ex.Y))
bounds2extent(bds) = Extent(X=bds[1], Y=bds[2])


function process_access_groups(output_folder, access_groups, boa_cube, qai_cube)
    # warp flags from cube to dop
    warp_flags = Dict(
                :s_srs=>convert(String, crs(boa_cube)),
                :t_srs=>convert(String, crs(first(access_groups).r[1])),
                )

    for (group_idx, group) in enumerate(access_groups)
        for tnr in unique(group.tnr)
            mkpath(joinpath(output_folder, tnr))
        end

        # get block xy
        x = first(first(group.xr))
        y = first(first(group.yr))
        boa_tile = boa_cube[y,x]
        qai_tile = qai_cube[y,x]

        block_groups = groupby(group, :blockrange)
        # process every time step independently
        for (time, boa, qai) in zip(string.(dims(boa_tile,Ti)), boa_tile, qai_tile)
            for blockgroup in block_groups
                try
                    dops = blockgroup.r
                    tnrs = blockgroup.tnr
                    # calculate the combined extent of all rasters on this tile / block combination
                    total_extent = combined_extent(dops, crs(boa))

                    local_qai_crop = qai[total_extent]
                    
                    dops_with_data = eltype(dops)[]
                    tnrs_with_data = eltype(tnrs)[]
                    qais_with_data = Raster[]
                    for (tnr, dop) in zip(tnrs, dops)
                        ex = reprojected_extent_with_margin(dop, crs(boa))
                        local_qai = local_qai_crop[ex]
                        # apply bitmask sets garbage values to 1 which equals the missingval and therefore can be filtered
                        mask = FCA.apply_bitmask(local_qai, cloudy_or_nodata)
                        has_data = FCA.contains_data(mask)
                        if has_data
                            push!(dops_with_data, dop)
                            push!(tnrs_with_data, tnr)
                            push!(qais_with_data, mask)
                            qai_out = warp(local_qai, warp_flags)[extent_with_margin(dop, 10)]
                            outname_qai = "$(tnr)_$(time[1:10])_QAI.tif"
                            write(joinpath(output_folder, tnr, outname_qai), qai_out; force=true)
                        end
                    end
                    
                    if length(dops_with_data) == 0
                        continue
                    end

                    data_extent = combined_extent(dops_with_data, crs(boa))
                    local_boa_crop = boa[data_extent]
                    
                    for (tnr, dop, mask) in zip(tnrs_with_data, dops_with_data, qais_with_data)
                        ex = reprojected_extent_with_margin(dop, crs(boa))
                        local_boa = local_boa_crop[ex]
                        @views for band in eachband(local_boa)
                            band[mask] .= band.missingval
                        end
                        boa_out = warp(local_boa, warp_flags)[extent_with_margin(dop, 10)]
                        outname_boa = "$(tnr)_$(time[1:10])_BOA.tif"
                        write(joinpath(output_folder, tnr, outname_boa), boa_out; force=true)
                    end
                catch err
                    @warn "Group with index $group_idx failed!"
                    @warn err
                end
            end
        end
    end
end


function process_overlapping_plots(boa_cube, qai_cube, rasters, output_folder)
    raster_crs = crs(first(rasters))
    warp_flags = Dict(
                :s_srs=>convert(String, crs(boa_cube)),
                :t_srs=>convert(String, raster_crs),
                )
    
    for r in rasters
        fname = basename(Rasters.filename(r))
        fname = split(fname, ".")[1]
        println(fname)

        try
            # get extent
            ((xmin, xmax),(ymin, ymax)) = Extents.extent(r, (X,Y))
            pts = [(xmin, ymin),(xmax, ymin),(xmax, ymax),(xmin, ymax)]
            # reproject from raster crs to cube crs
            bds = AG.reproject(pts, raster_crs, crs(boa_cube); order=:trad)
            # find bounding box and add margin
            xs = first.(bds)
            ys = last.(bds)
            xmin_, xmax_ = extrema(xs)
            ymin_, ymax_ = extrema(ys)
            s = step(dims(boa_cube,X))  # assume square pixels
            xmin_ -= 3*s
            xmax_ += 3*s
            ymin_ -= 3*s
            ymax_ += 3*s
            # preselect the qai data
            crop = qai_cube[X(xmin_..xmax_), Y(ymin_..ymax_)]

            if length(get_data(crop)) == 0
                @warn "Skipping $fname because it is not in the cube."
                continue
            end

            crop = read(crop)
            qai_masked = apply_bitmask(crop, cloudy_or_nodata)  # garbage is 1, usable data is 0, missingval should be 1
            qai_masked = extract_nonmissing(qai_masked)  # removes cutouts with only missingval
            
            # apply the selection to the BOA cube and load the data from disk
            # indexing with another cube produces a RasterSeries
            boa_series  = read.(boa_cube[qai_masked])
            qai_series  = read.(crop[qai_masked])
            mask_series = read.(seriesrepresentation(qai_masked))

            if length(boa_series) != length(mask_series)
                @warn "Extracted series have differing length for $fname."
            end
            
            # # set garbage values to nodata
            # # this should be packed into a function
            for (tile, mask_) in zip(boa_series, mask_series)
                for band in eachband(tile)
                    band[mask_] .= missingval(tile)
                end
            end
            
            boa_series = map(boa_series) do raster
                target_raster = warp(raster, warp_flags)
                target_raster[X(xmin-s..xmax+s), Y(ymin-s..ymax+s)]
            end
            
            qai_series = map(qai_series) do raster
                target_raster = warp(raster, warp_flags)
                target_raster[X(xmin-s..xmax+s), Y(ymin-s..ymax+s)]
            end
            
            mkpath(joinpath(output_folder, fname))
            
            for (i, file) in enumerate(boa_series)
                time = string(dims(boa_series, Ti)[i])
                outname_boa = fname * "_" * time[1:10] * "_BOA.tif"
                write(joinpath(output_folder, fname, outname_boa), file; force=true)
            end
            
            for (i, file) in enumerate(qai_series)
                time = string(dims(qai_series, Ti)[i])
                outname_qai = fname * "_" * time[1:10] * "_QAI.tif"
                write(joinpath(output_folder, fname, outname_qai), file; force=true)
            end

        catch err
            if isa(err, GDAL.GDALError)
                @warn "Couldn't read tile $fname due to a GDAL error while reading from the FORCE cube."
                continue
            else
                @warn "File $fname errored"
                @warn err
            end
        end
    end
end
#%%


function main(worker_id, world_size)
    boa_file = "boa_2023-10-03.jls191"
    qai_file = "qai_2023-10-03.jls191"
    output_folder = "/mnt/data/cutouts_2023"
    rasters  = deserialize("raster_metadata_2023.jls191");
    boa_cube = deserialize(boa_file);
    qai_cube = deserialize(qai_file);

    access_groups, overlapping_plots = get_access_groups(qai_cube, rasters);
    access_groups = groupby(access_groups, [:xr, :yr])

    access_group_chunk = chunks(1:length(access_groups), world_size)[worker_id][1]
    overlapping_chunk  = chunks(1:nrow(overlapping_plots), world_size)[worker_id][1]

    @time process_overlapping_plots(boa_cube, qai_cube, overlapping_plots.r[overlapping_chunk], output_folder)
    @time process_access_groups(output_folder, access_groups[access_group_chunk], boa_cube, qai_cube)
end

@time main(parse(Int, ARGS[1]), parse(Int, ARGS[2]))