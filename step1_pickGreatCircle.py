#!/usr/bin/env python3

# Use the proper idiom in the main module ...
# NOTE: See https://docs.python.org/3.12/library/multiprocessing.html#the-spawn-and-forkserver-start-methods
if __name__ == "__main__":
    # Import standard modules ...
    import argparse
    import gzip
    import json
    import math
    import os

    # Import special modules ...
    try:
        import cartopy
        cartopy.config.update(
            {
                "cache_dir" : os.path.expanduser("~/.local/share/cartopy_cache"),
            }
        )
    except:
        raise Exception("\"cartopy\" is not installed; run \"pip install --user Cartopy\"") from None
    try:
        import geojson
        geojson.geometry.Geometry.__init__.__defaults__ = (None, False, 12)     # NOTE: See https://github.com/jazzband/geojson/issues/135#issuecomment-596509669
    except:
        raise Exception("\"geojson\" is not installed; run \"pip install --user geojson\"") from None
    try:
        import matplotlib
        matplotlib.rcParams.update(
            {
                       "backend" : "Agg",                                       # NOTE: See https://matplotlib.org/stable/gallery/user_interfaces/canvasagg.html
                    "figure.dpi" : 300,
                "figure.figsize" : (9.6, 7.2),                                  # NOTE: See https://github.com/Guymer/misc/blob/main/README.md#matplotlib-figure-sizes
                     "font.size" : 8,
            }
        )
        import matplotlib.pyplot
    except:
        raise Exception("\"matplotlib\" is not installed; run \"pip install --user matplotlib\"") from None
    try:
        import numpy
    except:
        raise Exception("\"numpy\" is not installed; run \"pip install --user numpy\"") from None
    try:
        import shapely
        import shapely.geometry
        import shapely.ops
        import shapely.wkb
    except:
        raise Exception("\"shapely\" is not installed; run \"pip install --user Shapely\"") from None

    # Import my modules ...
    try:
        import pyguymer3
        import pyguymer3.geo
        import pyguymer3.image
    except:
        raise Exception("\"pyguymer3\" is not installed; you need to have the Python module from https://github.com/Guymer/PyGuymer3 located somewhere in your $PYTHONPATH") from None

    # **************************************************************************

    # Create argument parser and parse the arguments ...
    parser = argparse.ArgumentParser(
           allow_abbrev = False,
            description = "Show all the great circles from the starting location to it's antipode and choose a pair of them.",
                 epilog = "You cannot just draw the great circle to the antipode because there are an infinite number of possibilities: you can depart in any direction and you will end up in the same place eventually. Therefore, to get a fan of great circles going in all directions you must aim for a ring of points around the destination rather than the destination itself.",
        formatter_class = argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--altitude",
        default = 35.0e3,
           dest = "alt",
           help = "the altitude to fly above the chosen pair of great circles (in feet)",
           type = float,
    )
    parser.add_argument(
        "--debug",
        action = "store_true",
          help = "print debug messages",
    )
    parser.add_argument(
        "--maximum-segment-length",
        default = 5.0,
           dest = "maxSegLen",
           help = "the maximum length of any segment along a great circle (in nautical miles)",
           type = float,
    )
    parser.add_argument(
        "--number-of-great-circles",
        default = 9,
           dest = "nAngs",
           help = "the number of great circles to survey",
           type = int,
    )
    parser.add_argument(
        "--radius",
        default = 25.0,
           dest = "rad",
           help = "the radius around the antipode to stop at before actually arriving at the antipode (in nautical miles)",
           type = float,
    )
    parser.add_argument(
        "--starting-latitude",
        default = 50.5,
           dest = "startLat",
           help = "the latitude of the starting location",
           type = float,
    )
    parser.add_argument(
        "--starting-longitude",
        default = -1.0,
           dest = "startLon",
           help = "the longitude of the starting location",
           type = float,
    )
    parser.add_argument(
        "--timeout",
        default = 60.0,
           help = "the timeout for any requests/subprocess calls (in seconds)",
           type = float,
    )
    parser.add_argument(
        "--tolerance",
        default = 1.0e-10,
           dest = "tol",
           help = "the Euclidean distance that defines two points as being the same (in degrees)",
           type = float,
    )
    args = parser.parse_args()

    # **************************************************************************

    # Make output directory if it is missing ...
    if not os.path.exists("output"):
        os.mkdir("output")

    # Save arguments ...
    with open("output/args.json", "wt", encoding = "utf-8") as fObj:
        json.dump(
            vars(args),
            fObj,
            ensure_ascii = False,
                  indent = 4,
               sort_keys = True,
       )

    # Define starting location and antipode ...
    startPnt = shapely.geometry.point.Point(args.startLon, args.startLat)
    antipPnt = shapely.geometry.point.Point(args.startLon + 180.0, -args.startLat)

    # Convert user inputs to useful units and then figure out how far away the
    # horizon is from the specified altitude ...
    args.alt *= 12.0 * 0.0254                                                   # [m]
    args.rad *= 1852.0                                                          # [m]
    args.maxSegLen *= 1852.0                                                    # [m]
    pathWidth = pyguymer3.consts.RADIUS_OF_EARTH * math.acos(pyguymer3.consts.RADIUS_OF_EARTH / (pyguymer3.consts.RADIUS_OF_EARTH + args.alt))  # [m]

    print(f"The horizon is {0.001 * pathWidth:,.1f} km away (along the surface of the Earth) when viewed from an altitude of {0.001 * args.alt:,.1f} km.")

    # Find the buffer around the antipode ...
    antipPly = pyguymer3.geo.buffer(
        antipPnt,
        args.rad,
        debug = args.debug,
         fill = -1.0,
         nang = args.nAngs,
         simp = -1.0,
          tol = args.tol,
    )

    # Initialize array ...
    dists = numpy.zeros(args.nAngs, dtype = numpy.float64)                      # [m]

    # Loop over angles around the antipode ...
    for iAng in range(args.nAngs):
        print(f"Calculating the distance to point {iAng + 1:d}/{args.nAngs:d} around the antipode ...")

        # Create short-hands ...
        antipLon, antipLat = antipPly.exterior.coords[iAng]                     # [°], [°]

        # Calculate the distance to the point around the antipode ...
        dists[iAng], _, _ = pyguymer3.geo.calc_dist_between_two_locs(
            args.startLon,
            args.startLat,
            antipLon,
            antipLat,
        )                                                                       # [m]

    # Create short-hands ...
    minDist = dists.min()                                                       # [m]
    maxDist = dists.max()                                                       # [m]

    print(f"The minimum distance to within {0.001 * args.rad:,.1f} km of the antipode is {0.001 * minDist:,.1f} km.")
    print(f"The maximum distance to within {0.001 * args.rad:,.1f} km of the antipode is {0.001 * maxDist:,.1f} km.")

    # Pick the two longest great circles ...
    # NOTE: The Earth is fatter over the equator than over the poles. Therefore,
    #       great circles which spend more time near the equator (i.e., great
    #       circles with more horizontal paths through the equatorial region)
    #       are longer than great circles which spend more time near the poles
    #       (which will intersect the equator more vertically).
    key1, key2 = sorted(dists.argsort()[::-1][:2])                              # [#]

    print(f"The two chosen great circles are for angles {key1 + 1:d}/{args.nAngs:d} and {key2 + 1:d}/{args.nAngs:d}.")

    # **************************************************************************

    # Create figure ...
    fg = matplotlib.pyplot.figure()

    # Create axes ...
    axL = pyguymer3.geo.add_axis(
        fg,
        coastlines_resolution = "i",
                        debug = args.debug,
                         dist = 1.5 * args.rad,
                        index = 1,
                          lat = args.startLat,
                          lon = args.startLon,
                        ncols = 2,
                        nrows = 2,
                          tol = args.tol,
    )
    axR = pyguymer3.geo.add_axis(
        fg,
        coastlines_resolution = "i",
                        debug = args.debug,
                         dist = 1.5 * args.rad,
                        index = 2,
                          lat = -args.startLat,
                          lon = args.startLon + 180.0,
                        ncols = 2,
                        nrows = 2,
                          tol = args.tol,
    )
    axB = pyguymer3.geo.add_axis(
        fg,
        coastlines_resolution = "c",
                        debug = args.debug,
                        index = (3, 4),
                        ncols = 2,
                        nrows = 2,
                          tol = args.tol,
    )

    # Configure axes ...
    pyguymer3.geo.add_map_background(
        axL,
             debug = args.debug,
              name = "shaded-relief",
        resolution = "large8192px",
    )
    pyguymer3.geo.add_map_background(
        axR,
             debug = args.debug,
              name = "shaded-relief",
        resolution = "large8192px",
    )
    pyguymer3.geo.add_map_background(
        axB,
             debug = args.debug,
              name = "shaded-relief",
        resolution = "large8192px",
    )

    # Initialize list ...
    lines = []

    # Loop over angles around the antipode ...
    for iAng in range(args.nAngs):
        print(f"Calculating great circle for angle {iAng + 1:d}/{args.nAngs:d} ...")

        # Create short-hands ...
        antipLon, antipLat = antipPly.exterior.coords[iAng]                     # [°], [°]

        # Find the great circle to the point around the antipode ...
        circle = pyguymer3.geo.great_circle(
            args.startLon,
            args.startLat,
            antipLon,
            antipLat,
              debug = args.debug,
            maxdist = args.maxSegLen,
        )

        # Check if this great circle is a chosen one ...
        if iAng in [key1, key2]:
            # Add the lines of this great circle to the list ...
            lines += pyguymer3.geo.extract_lines(circle)

        # Draw the great circle and clean up ...
        axL.add_geometries(
            pyguymer3.geo.extract_lines(circle),
            cartopy.crs.PlateCarree(),
            edgecolor = matplotlib.pyplot.cm.jet((dists[iAng] - minDist) / (maxDist - minDist)),
            facecolor = "none",
            linewidth = 1.0,
               zorder = 5.1,
        )
        axR.add_geometries(
            pyguymer3.geo.extract_lines(circle),
            cartopy.crs.PlateCarree(),
            edgecolor = matplotlib.pyplot.cm.jet((dists[iAng] - minDist) / (maxDist - minDist)),
            facecolor = "none",
            linewidth = 1.0,
               zorder = 5.1,
        )
        axB.add_geometries(
            pyguymer3.geo.extract_lines(circle),
            cartopy.crs.PlateCarree(),
            edgecolor = matplotlib.pyplot.cm.jet((dists[iAng] - minDist) / (maxDist - minDist)),
            facecolor = "none",
            linewidth = 1.0,
               zorder = 5.1,
        )
        del circle

    # Clean up ...
    del antipPly

    # **************************************************************************

    print("Calculating path (as a line) ...")

    # Create a sorted array of all of the coordinates in the great circles
    # joined together from the anti-meridian at -180.0° to the anti-meridian at
    # +180.0° (remembering to add the antipode itself so that there isn't a gap
    # where all the great circles meet) ...
    tmpCoords = []                                                              # [°], [°]
    for line in lines:
        for coord in line.coords:
            if coord in tmpCoords:
                print(f"INFO: ({coord[0]:e}°,{coord[1]:e}°) is already in \"tmpCoords\".")
                continue
            tmpCoords.append(coord)                                             # [°], [°]
    tmpCoords.append((args.startLon + 180.0, -args.startLat))                   # [°], [°]
    tmpCoords = numpy.array(tmpCoords)                                          # [°]
    coords = numpy.zeros(tmpCoords.shape, dtype = numpy.float64)                # [°]
    for i, key in enumerate(tmpCoords[:, 0].argsort()):
        coords[i, :] = tmpCoords[key, :]                                        # [°]
    del tmpCoords

    print("Saving path (as a line) ...")

    # Save line and clean up ...
    with open("output/pathLine.csv", "wt", encoding = "utf-8") as fObj:
        fObj.write("longitude [°],latitude [°]\n")
        for i in range(coords.shape[0]):
            fObj.write(f"{coords[i, 0]:e},{coords[i, 1]:e}\n")
    del coords

    print("Calculating path (as a band) ...")

    # Initialize list ...
    polys = []

    # Loop over lines ...
    for line in lines:
        # Add the buffer of this line to the list ...
        polys += pyguymer3.geo.extract_polys(
            pyguymer3.geo.buffer(
                line,
                pathWidth,
                debug = args.debug,
                 fill = -1.0,
                 nang = args.nAngs,
                 simp = -1.0,
                  tol = args.tol,
            )
        )

    # Clean up ...
    del lines

    # Convert list of Polygons to a (unified) [Multi]Polygon and clean up ...
    path = shapely.ops.unary_union(polys).simplify(args.tol)
    del polys

    print("Saving path (as a band) ...")

    # Save [Multi]Polygon ...
    with gzip.open("output/pathBand.wkb.gz", mode = "wb", compresslevel = 9) as gzObj:
        gzObj.write(shapely.wkb.dumps(path))

    # Save [Multi]Polygon ...
    with open("output/pathBand.geojson", "wt", encoding = "utf-8") as fObj:
        geojson.dump(
            path,
            fObj,
            ensure_ascii = False,
                  indent = 4,
               sort_keys = True,
        )

    # Draw the path and clean up ...
    axB.add_geometries(
        pyguymer3.geo.extract_polys(path),
        cartopy.crs.PlateCarree(),
        edgecolor = "none",
        facecolor = (1.0, 0.0, 0.0, 0.5),
        linewidth = 1.0,
           zorder = 5.0,
    )
    del path

    # **************************************************************************

    # Configure figure ...
    fg.suptitle(f"Great circles from ({args.startLat:+.1f}°,{args.startLon:+.1f}°) to a ring {0.001 * args.rad:,.1f} km around its antipode")
    fg.tight_layout()

    print("Saving figure ...")

    # Save figure ...
    fg.savefig("output/greatCircles.png")
    matplotlib.pyplot.close(fg)

    # Optimise PNG ...
    pyguymer3.image.optimize_image(
        "output/greatCircles.png",
          debug = args.debug,
          strip = True,
        timeout = args.timeout,
    )
