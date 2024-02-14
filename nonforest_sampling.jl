# load tree cover map
# reproject crs
# load and calculate the extents where points should be sampled
# sample random points within extents -> accept if in forest
# constrain to area of germany if necessary
# save the points in some way
# add to list with points to extract

using GDAL
using Base.Iterators
using ArchGDAL
import ArchGDAL as AG
using Rasters


function reprojected_extent(r::Raster, target_crs)
    ((xmin, xmax),(ymin, ymax)) = Rasters.Extents.extent(r, (X,Y))
    pts = [(xmin, ymin),(xmax, ymin),(xmax, ymax),(xmin, ymax)]
    # reproject from raster crs to target crs
    bds = AG.reproject(pts, crs(r), target_crs; order=:trad)
    xs = first.(bds)
    ys = last.(bds)
    xmin_, xmax_ = extrema(xs)
    ymin_, ymax_ = extrema(ys)
    return Rasters.Extents.Extent((X=(xmin_, xmax_), Y=(ymin_, ymax_)))
end


function random_points(ex, n, min_dist=14.2)
    out = Vector{Tuple{Float64, Float64}}(undef, n)
    (xmin, xmax), (ymin, ymax) = ex.X, ex.Y
    x = rand() * (xmax - xmin) + xmin
    y = rand() * (ymax - ymin) + ymin
    out[1] = (x, y)
    
    i = 2
    while i <= n
        x = rand() * (xmax - xmin) + xmin
        y = rand() * (ymax - ymin) + ymin
        distances = (((x-x2)^2 + (y-y2)^2) > min_dist^2 for (x2, y2) in view(out, 1:i-1))
        if all(distances)
            out[i] = (x, y)
            i += 1
        end
    end
    return out
end


function write_result_sqlite(filename, result)
    AG.create(filename; driver=AG.getdriver("SQLite")) do ds
        AG.createlayer(name = "geom",
                       dataset = ds,
                       geom = AG.wkbPoint,
                       spatialref = AG.importEPSG(25832)) do layer
            AG.addfielddefn!(layer, "tree_id", AG.OFTInteger)
            AG.addfielddefn!(layer, "plot_id", AG.OFTInteger)
            id_idx = AG.findfieldindex(AG.layerdefn(layer), "tree_id")
            plot_id_idx = AG.findfieldindex(AG.layerdefn(layer), "plot_id")
            i = 0
            for chunk in Iterators.partition(result, 20000)
                GDAL.ogr_l_starttransaction(layer)
                for r in chunk
                    id, pts = r
                    for p in pts
                        AG.addfeature(layer) do feature
                            AG.setfid!(feature, i)
                            AG.setfield!(feature, id_idx, -i-1)
                            AG.setfield!(feature, plot_id_idx, id)
                            AG.setgeom!(feature, 0, AG.createpoint(p...))
                            return nothing
                        end
                        i += 1 
                    end
                end
                GDAL.ogr_l_committransaction(layer)
            end
        end
    end
    return nothing
end


function main(tcm_file, dop_dir, max_points_per_plot=6, samples_per_plot=15)
    tcm = Raster(tcm_file; lazy=true)
    dops = [Raster(f; lazy=true) for f in readdir(dop_dir; join=true)]
    plot_ids = parse.(Int, first.(splitext.(basename.((Rasters.filename.(dops))))))
    out = Tuple{Int, Vector{Tuple{Float64, Float64}}}[]
    sizehint!(out, 27000*5)

    println("start generating points")
    @time for (id, dop) in zip(plot_ids, dops)
        rep_ex = reprojected_extent(dop, crs(tcm))
        no_trees = tcm[rep_ex] .< 0.1
        pts = random_points(rep_ex, samples_per_plot)
        pts = filter(pts) do point
            x,y = point
            return all(no_trees[X(x-20..x+20), Y(y-20..y+20)])
        end
        pts = Tuple.(AG.reproject(pts, crs(tcm), EPSG(25832); order=:trad))
        if length(pts) > 0
            push!(out, (id, pts[1:min(length(pts), max_points_per_plot)]))
        end
    end
    return out
end

#%%
tcm_file = "/data_hdd/treecover_germany_2018.tif"
dop_dir = "/data_hdd/bkg/clip_2023/"
@time res = main(tcm_file, dop_dir)

#%%
write_result_sqlite("nonforest_pixels.sqlite", res)
