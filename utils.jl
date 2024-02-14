using Dates
using SQLite
import ArchGDAL as AG
import GDAL
using DataFrames

printlayers(ds) =
    for i in 0:(AG.nlayer(ds) - 1)
        println("$i: $(AG.getname(AG.getlayer(ds, i)))")
    end

function printfields(layer)
    for i in 0:(AG.nfield(layer) - 1)
        println("$i: $(AG.getname(AG.getfielddefn(AG.layerdefn(layer), i)))")
    end
end

function printds(ds)
    for i in 0:(AG.nlayer(ds) - 1)
        layer = AG.getlayer(ds, i)
        println("$i: $(AG.getname(layer))")

        for j in 0:(AG.nfield(layer) - 1)
            println("\t$j: $(AG.getname(AG.getfielddefn(AG.layerdefn(layer), j)))")
        end
    end
end

parseWKB(wkb_hex) = ismissing(wkb_hex) ? missing : AG.fromWKB(hex2bytes(wkb_hex))

function clearspatialfilter(layer::AG.IFeatureLayer)
    return GDAL.ogr_l_setspatialfilter(layer.ptr, C_NULL)
end

function executesql(f,
        dataset::AG.AbstractDataset,
        query::AbstractString;
        dialect = "",
        spatialfilter = AG.Geometry(C_NULL))
    ret = AG.unsafe_executesql(dataset, query; dialect, spatialfilter)
    res = f(ret)
    AG.releaseresultset(dataset, ret)
    return res
end

function executesql(dataset::AG.AbstractDataset,
        query::AbstractString;
        dialect = "",
        spatialfilter = AG.Geometry(C_NULL))
    return executesql(DataFrame, dataset, query; dialect = dialect,
        spatialfilter = spatialfilter)
end

flatten(a) = map(x -> getindex.(a, x), 1:length(first(a)))

rescale(x, xmin, xmax, ymin, ymax) = (ymax - ymin) / (xmax - xmin) * (x - xmin) + ymin

function rescale(p::CartesianIndex{2}, ncols, nrows, xr_to, yr_to)
    y, x = Tuple(p)
    x2 = rescale(x, 1, ncols, extrema(xr_to)...)
    y2 = rescale(y, 1, nrows, extrema(yr_to)...)
    return y2, x2
end

skip_iter(v, i) = (v[j] for j in firstindex(v):lastindex(v) if !(j in i))

function serialize_df(fname, df)
    local_df = copy(df)
    local_df.tree_position = GI.coordinates.(local_df.tree_position)
    local_df.disk = GI.coordinates.(local_df.disk)
    local_df.plot_center = GI.coordinates.(local_df.plot_center)
    return serialize(fname, local_df)
end

function deserialize_df(fname)
    df = deserialize(fname)
    df.tree_position = AG.createpoint.(df.tree_position)
    df.disk = AG.createpolygon.(df.disk)
    df.plot_center = AG.createpoint.(df.plot_center)
    return df
end

function envelopetopolygon(envelope)
    xmin = envelope.MinX
    ymin = envelope.MinY
    xmax = envelope.MaxX
    ymax = envelope.MaxY
    return AG.createpolygon([(xmin, ymin), (xmax, ymin), (xmax, ymax), (xmin, ymax),
        (xmin, ymin)])
end

function writepolygonstodisk(filepath, polygons)
    # gdal_polygons = convert.(AG.IGeometry{AG.wkbPolygon}, polygons)
    layername = basename(splitext(filepath)[1])
    AG.create(AG.getdriver("Memory")) do ds
        AG.createlayer(; name = layername, dataset = ds, geom = AG.wkbUnknown) do layer
            for polygon in polygons
                AG.createfeature(layer) do feature
                    return AG.setgeom!(feature, polygon)
                end
            end
        end
        return AG.write(ds,
            filepath;
            driver = AG.getdriver("SQLite"),
            use_gdal_copy = false)
    end
end

function read_force_series(path, tnr)
    boas = filter(contains("BOA"), readdir(joinpath(path, string(tnr)); join = true))
    qais = filter(contains("QAI"), readdir(joinpath(path, string(tnr)); join = true))
    times = DateTime.(getindex.(split.(basename.(boas), '_'), 2))
    boa_series = RasterSeries(boas, Ti(times); duplicate_first = true, lazy = true)
    qai_series = RasterSeries(qais, Ti(times); duplicate_first = true, lazy = true)
    return boa_series, qai_series
end

vec2blob(v) = String(reinterpret(UInt8, v))

function prepare_df_dtypes_for_export!(df)
    if eltype(df.boa) != String
        df.boa = map(vec2blob, df.boa)
    end
    if !("doy" in names(df))
        df.doy = Int32.(dayofyear.(df.time))
    end
    if eltype(df.time) <: DateTime
        # remember to update this before 2038 ;-)
        df.time = Int32.(datetime2unix.(df.time))
    end
    df.tree_id = Int32.(df.tree_id)
    df.qai = Int32.(df.qai)
    return df
end

function write_sqlite(fpath, df)
    dbname = "data" #first(splitext(basename(fpath)))
    db = SQLite.DB(fpath)
    SQLite.load!(df, db, dbname)
    close(db)
    return nothing
end

function import_df_dtypes!(df)
    df.boa = reinterpret.(Int16, codeunits.(df.boa))
    df.time = Dates.unix2datetime.(df.time)
    GC.gc()
    return df
end

"""
    valid_polygon_indices(polygons::T, cutoff, distance) where T

Checks which tree crown projection circles are valid and returns their indices.
The check is performed as follows:
- For each tree, the distance to all other trees, as well as the overlap of 
the projection circle with the union of all other projection circles is calculated.
- Then it is checked, whether the current tree is the biggest (based on polygon area) 
within `distance`.
- If the tree is not the biggest, has other trees within `distance` and other trees 
    cover more than overlap_cutoff of its projection circle, the tree is omitted.
- Otherwise it is valid
"""
function valid_polygon_indices(polygons::T, overlap_cutoff, distance) where {T}
    N = length(polygons)
    N == 1 && return [firstindex(polygons)]
    valid_polygons = Int[]
    for i in 1:N
        polygon = polygons[i]
        center = AG.centroid(polygon)
        area = AG.geomarea(polygon)
        other_polygons = skip_iter(polygons, i)
        other_areas = AG.geomarea.(other_polygons)

        other_centroids = AG.centroid.(other_polygons)
        distances = [AG.distance(center, oc) for oc in other_centroids]
        close_trees = distances .< distance
        is_biggest = all(<(area), view(other_areas, close_trees))

        if is_biggest || (sum(close_trees) == 0)
            push!(valid_polygons, i)
            continue
        end

        other_union = foldl(AG.union,
            other_polygons)::Union{AG.IGeometry{AG.wkbPolygon},
            AG.IGeometry{AG.wkbMultiPolygon}}
        intersection_area = AG.geomarea(AG.intersection(polygon, other_union))
        covered_fraction = intersection_area / area

        if covered_fraction < overlap_cutoff
            push!(valid_polygons, i)
        end
    end
    return valid_polygons
end