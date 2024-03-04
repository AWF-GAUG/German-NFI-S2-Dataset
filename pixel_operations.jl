using Rasters
import ArchGDAL as AG

"""
    underlying_pixelgrid(polygon, xmin_global, ymin_global, xstep, ystep)

Creates a point grid underlying the given polygon. Arguments are the polygon 
and the global x,y origin of an underlying image as well as the resolution. 
The resulting grid coordinates are the origins of all pixels covered by the 
polygon. Only tested for ETRS25832 (Germany).
"""
function underlying_pixelgrid(polygon, xmin_global, ymin_global, xstep, ystep)
    # we must work in the coordinate system of the Rasters package
    # pixel origin is south-west corner in ETRS25832
    xstep > 0 || error("pixel grid generation only works for positive x step")
    ystep < 0 || error("pixel grid generation only works for negative y step")
    envelope = AG.envelope(polygon)
    xmin = envelope.MinX - xmin_global
    xmax = envelope.MaxX - xmin_global
    ymin = envelope.MinY - ymin_global
    ymax = envelope.MaxY - ymin_global
    xpxmin = floor(xmin / xstep) * xstep + xmin_global
    xpxmax = floor(xmax / xstep) * xstep + xmin_global
    ypxmin = ceil(ymin / ystep) * ystep + ymin_global  # use ceil because arg is negative
    ypxmax = ceil(ymax / ystep) * ystep + ymin_global
    xrange = xpxmin:abs(xstep):(xpxmax + 1e-6)  # add a little to avoid floating point stuff
    yrange = ypxmin:abs(ystep):(ypxmax + 1e-6)
    return tuple.(xrange', yrange)
end

"""
    coord_to_polygon(coord, xstep, ystep)

Takes a (xmin, ymin) coordinate and the resolution in x and y direction
and returns a polygon covering the respective pixel.
"""
function coord_to_polygon(coord, xstep, ystep)
    xmin, ymin = coord
    AG.createpolygon([(xmin, ymin),
        (xmin + xstep, ymin),
        (xmin + xstep, ymin - ystep),
        (xmin, ymin - ystep),
        (xmin, ymin)])
end

"""
    mean_pixel_value(polygon, pixel_coords, pixel_polygons, raster)

Computes the area-weighted average of all pixels intersected by a polygon.

Args:
- polygon: The polygon
- pixel_coords: The coordinates of the underlying pixel grid
- pixel_polygons: The pixel polygons of the underlying image
- raster: The source raster from which the values are extracted
"""
function weighted_mean_pixel_value(polygon, pixel_coords, pixel_polygons, raster)
    pa = AG.geomarea(polygon)
    res = zeros(Float32, 10)
    for (coord, pixel) in zip(pixel_coords, pixel_polygons)
        x, y = coord
        pixel_value = @view raster[X(Near(x)), Y(Near(y))]
        
        # return nodata if there's a nodata pixel in the area
        if pixel_value[1] == -9999
            res[:] .= -9999
            return Int16.(res)
        end
        
        pixel_weight = AG.geomarea(AG.intersection(polygon, pixel)) / pa
        res += pixel_value * pixel_weight
    end
    return round.(Int16, res)
end