#!/usr/bin/env python3

# Use the proper idiom in the main module ...
# NOTE: See https://docs.python.org/3.11/library/multiprocessing.html#the-spawn-and-forkserver-start-methods
if __name__ == "__main__":
    # Import standard modules ...
    import argparse
    import gzip
    import os
    import shutil

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
        import shapely.wkb
    except:
        raise Exception("\"shapely\" is not installed; run \"pip install --user Shapely\"") from None

    # Import my modules ...
    try:
        import pyguymer3
        import pyguymer3.geo
        import pyguymer3.image
        import pyguymer3.media
    except:
        raise Exception("\"pyguymer3\" is not installed; you need to have the Python module from https://github.com/Guymer/PyGuymer3 located somewhere in your $PYTHONPATH") from None

    # **************************************************************************

    # Create argument parser and parse the arguments ...
    parser = argparse.ArgumentParser(
           allow_abbrev = False,
            description = "Animate the journey along the chosen great circle and every country that it intersects.",
        formatter_class = argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--debug",
        action = "store_true",
          help = "print debug messages",
    )
    parser.add_argument(
        "--number-of-frames",
        default = 120,
           dest = "nFrames",
           help = "the number of frames in the animation",
           type = int,
    )
    parser.add_argument(
        "--radius-of-arrows",
        default = 0.1 * pyguymer3.consts.CIRCUMFERENCE_OF_EARTH,
           dest = "radOfArrows",
           help = "the Geodetic distance from the centre of the map to the limit of bothering to draw the arrows (in metres)",
           type = float,
    )
    parser.add_argument(
        "--radius-of-labels",
        default = 0.2 * pyguymer3.consts.CIRCUMFERENCE_OF_EARTH,
           dest = "radOfLabels",
           help = "the Geodetic distance from the centre of the map to the labels (in metres)",
           type = float,
    )
    parser.add_argument(
        "--timeout",
        default = 60.0,
           help = "the timeout for any requests/subprocess calls",
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

    # Load the path and process it ...
    pathLons, pathLats = numpy.loadtxt(
        "output/pathLine.csv",
        delimiter = ",",
         skiprows = 1,
           unpack = True,
    )                                                                           # [°], [°]
    pathSegDists = numpy.zeros(pathLons.size - 1, dtype = numpy.float64)        # [m]
    for i in range(pathSegDists.size):
        pathSegDists[i], _, _ = pyguymer3.geo.calc_dist_between_two_locs(
            pathLons[i],
            pathLats[i],
            pathLons[i + 1],
            pathLats[i + 1],
        )                                                                       # [m]
    pathCumDists = numpy.zeros(pathLons.size, dtype = numpy.float64)            # [m]
    for i in range(1, pathCumDists.size):
        pathCumDists[i] = pathSegDists[:i].sum()                                # [m]

    # Load the path ...
    with gzip.open("output/pathBand.wkb.gz", mode = "rb") as gzObj:
        pathPolys = pyguymer3.geo.extract_polys(shapely.wkb.loads(gzObj.read()))

    # **************************************************************************

    # Initialize lists ...
    neNames = []
    intersectPolys = []

    # Find file containing all the country shapes ...
    sfile = cartopy.io.shapereader.natural_earth(
          category = "cultural",
              name = "admin_0_countries",
        resolution = "10m",
    )

    # Loop over records ...
    for record in cartopy.io.shapereader.Reader(sfile).records():
        # Create short-hand and merge territories into their parent country ...
        neName = pyguymer3.geo.getRecordAttribute(record, "NAME")
        if neName in ["Bajo Nuevo Bank", "Serranilla Bank",]:
            neName = "Colombia"
        if neName in ["U.S. Minor Outlying Is.", "USNB Guantanamo Bay",]:
            neName = "United States of America"

        # Loop over record Polygons ...
        for recordPoly in pyguymer3.geo.extract_polys(record.geometry):
            # Loop over path Polygons ...
            for pathPoly in pathPolys:
                # Calculate the intersection and skip if it is empty ...
                intersectMultiPoly = recordPoly.intersection(pathPoly)
                if intersectMultiPoly.is_empty:
                    continue

                # Loop over intersection Polygons ...
                for intersectPoly in pyguymer3.geo.extract_polys(intersectMultiPoly):
                    # Add intersection to lists ...
                    neNames.append(neName)
                    intersectPolys.append(intersectPoly)

    print(f"There are {len(intersectPolys):,d} intersection Polygons.")

    # Create short-hand ...
    uniqueNeNames = sorted(set(neNames))

    print("The following countries will be visited:")
    for uniqueNeName in uniqueNeNames:
        print(f"  * {uniqueNeName}")

    # **************************************************************************

    # Initialize lists ...
    midLats = []                                                                # [°]
    midLons = []                                                                # [°]

    # Loop over intersections ...
    for neName, intersectPoly in zip(neNames, intersectPolys):
        # Convert the CoordinateSequence of the exterior LinearRing of the
        # intersection into a NumPy array ...
        points = numpy.array(intersectPoly.exterior.coords)                     # [°]

        # Find the longitude of the middle of the points and clean up ...
        midLon, midLat, _ = pyguymer3.geo.find_middle_of_locs(
            points[:, 0],
            points[:, 1],
            debug = args.debug,
              pad = -1.0,
        )                                                                       # [°], [°]
        del points

        # Append values to lists ...
        midLats.append(midLat)                                                  # [°]
        midLons.append(midLon)                                                  # [°]

    # **************************************************************************

    # Initialize dictionary ...
    minMidLons = {}

    # Loop over the unique country names ...
    for uniqueNeName in uniqueNeNames:
        # Loop over intersections ...
        for neName, intersectPoly, midLon in zip(neNames, intersectPolys, midLons):
            # Skip this intersection if it is not from the correct country ...
            if neName != uniqueNeName:
                continue

            # Update dictionary ...
            if uniqueNeName not in minMidLons:
                minMidLons[uniqueNeName] = 1000.0                               # [°]
            minMidLons[uniqueNeName] = min(midLon, minMidLons[uniqueNeName])    # [°]

    # Find the order to loop through the country names to go from West to East ...
    sortedNeNames = sorted(minMidLons, key = minMidLons.get)

    print("Going from West to East these are:")
    for sortedNeName in sortedNeNames:
        print(f"  * {sortedNeName}")

    # **************************************************************************

    # Make output folder ...
    if not os.path.exists("output/frames"):
        os.mkdir("output/frames")

    # Create arrays to describe the journey ...
    fracDists = numpy.linspace(0.0, 1.0, endpoint = False, num = args.nFrames)[1:]
    dists = fracDists * pathSegDists.sum()                                      # [m]
    pngs = [f"output/frames/{dist:08.0f}.png" for dist in dists]

    # Loop over frames ...
    for dist, png in zip(dists, pngs):
        # Skip this frame if it already exists ...
        if os.path.exists(png):
            continue

        print(f"Making \"{png}\" ...")

        # Interpolate the latitude and longitude to find where this figure
        # should be centered ...
        figLat = numpy.interp(
             x = dist,
            xp = pathCumDists,
            fp = pathLats,
        )                                                                       # [°]
        figLon = numpy.interp(
             x = dist,
            xp = pathCumDists,
            fp = pathLons,
        )                                                                       # [°]

        # Create figure ...
        fg = matplotlib.pyplot.figure(figsize = (12.8, 12.8))

        # Create axis ...
        ax = pyguymer3.geo.add_axis(
            fg,
                            debug = args.debug,
            coastlines_resolution = "c",
                              lat = figLat,
                              lon = figLon,
                              tol = args.tol,
        )

        # Configure axis ...
        pyguymer3.geo.add_map_background(
            ax,
                 debug = args.debug,
                  name = "shaded-relief",
            resolution = "large8192px",
        )

        # Draw the path ...
        ax.add_geometries(
            pathPolys,
            cartopy.crs.PlateCarree(),
            edgecolor = "none",
            facecolor = (1.0, 0.0, 0.0, 0.25),
            linewidth = 1.0,
               zorder = 5.0,
        )

        # **********************************************************************

        # Loop over countries ...
        for i, sortedNeName in enumerate(sortedNeNames):
            # Label the intersection ...
            bear = 360.0 * float(i) / float(len(sortedNeNames))                 # [°]
            labelLon, labelLat, _ = pyguymer3.geo.calc_loc_from_loc_and_bearing_and_dist(
                figLon,
                figLat,
                bear,
                args.radOfLabels,
            )                                                                   # [°], [°]
            pyguymer3.geo.add_annotation(
                ax,
                labelLon,
                labelLat,
                sortedNeName,
                  bbox = {
                    "edgecolor" : "black",
                    "facecolor" : "white",
                    "linewidth" : 1.0,
                },
                 debug = args.debug,
                zorder = 5.3,
            )

            # Loop over intersections ...
            for neName, intersectPoly, midLat, midLon in zip(neNames, intersectPolys, midLats, midLons):
                # Skip this intersection if it is not from the correct country ...
                if neName != sortedNeName:
                    continue

                # Skip this intersection if it is too far away ...
                try:
                    if pyguymer3.geo.calc_dist_between_two_locs(
                        figLon,
                        figLat,
                        midLon,
                        midLat,
                    )[0] > args.radOfArrows:
                        continue
                except:
                    print(f"WARNING: Had to skip an intersection when calculating the distance between ({figLat:+.3f}°,{figLon:+.3f}°) and ({midLat:+.3f}°,{midLon:+.3f}°).")
                    continue

                # Draw the intersection ...
                ax.add_geometries(
                    [intersectPoly],
                    cartopy.crs.PlateCarree(),
                    edgecolor = "none",
                    facecolor = f"C{i % 10:d}",
                    linewidth = 0.0,
                       zorder = 5.1,
                )

                # Annotate the intersection ...
                pyguymer3.geo.add_annotation(
                    ax,
                    midLon,
                    midLat,
                    "",
                    arrowprops = {
                             "alpha" : 0.5,
                        "arrowstyle" : "-",
                         "edgecolor" : f"C{i % 10:d}",
                         "facecolor" : "none",
                         "linewidth" : 1.0,
                            "zorder" : 5.2,
                    },
                         debug = args.debug,
                        txtLat = labelLat,
                        txtLon = labelLon,
                )

        # **********************************************************************

        # Configure figure ...
        fg.suptitle(f"({figLat:+.3f}°,{figLon:+.3f}°) is {0.001 * dist:,.1f} km along the route")
        fg.tight_layout()

        # Save figure ...
        fg.savefig(png)
        matplotlib.pyplot.close(fg)

        # Optimise PNG ...
        pyguymer3.image.optimize_image(
            png,
              debug = args.debug,
              strip = True,
            timeout = args.timeout,
        )

    # **************************************************************************

    print("Making \"output/pathBand.mp4\" ...")

    # Save 60fps MP4 ...
    vname = pyguymer3.media.images2mp4(
        pngs,
        debug = args.debug,
          fps = 60.0,
    )
    shutil.move(vname, "output/pathBand.mp4")

    # Set heights ...
    # NOTE: By inspection, the PNG frames are 3,840 px wide.
    heights = [512, 1024, 2048]                                                 # [px]

    # Loop over heights ...
    for height in heights:
        print(f"Making \"output/pathBand{height:04d}px.mp4\" ...")

        # Save 60fps MP4 ...
        vname = pyguymer3.media.images2mp4(
            pngs,
                   debug = args.debug,
                     fps = 60.0,
             screenWidth = height,
            screenHeight = height,
        )
        shutil.move(vname, f"output/pathBand{height:04d}px.mp4")
