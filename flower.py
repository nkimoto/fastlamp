#!/usr/bin/env python3

"""
Copyright (c) 2013, LAMP development team
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * Neither the name of the LAMP development team nor the
      names of its contributors may be used to endorse or promote products
      derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL LAMP DEVELOPMENT TEAM BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""

import sys
import math
import argparse
import flower.flower_svg as svg
from flower.flower_readfile import readResult, combiName, motifName, motifNgenes, motifApvalue, combiApvalue, combiRank

__author__ = "Takayuki ITOH"


# definition of main function
def main():
    # add command line options
    parser = argparse.ArgumentParser(description="Generate flower-like visualizations from LAMP results.")
    parser.add_argument("result_filename", help="LAMP result file")
    parser.add_argument("target_filename", nargs='?', default='', help="Target (csv) filename (optional)")
    parser.add_argument("value_filename", nargs='?', default='', help="Value (tab) filename (optional)")
    parser.add_argument("--ewidth", dest="EWIDTH", type=int, default=120)
    parser.add_argument("--eheight", dest="EHEIGHT", type=int, default=50)
    parser.add_argument("--minwidth", dest="MINWIDTH", type=int, default=120)
    parser.add_argument("--minheight", dest="MINHEIGHT", type=int, default=40)
    parser.add_argument("--shiftx", dest="SHIFTX", type=int, default=400)
    parser.add_argument("--shifty", dest="SHIFTY", type=int, default=400)
    parser.add_argument("--petalposcoef", dest="PETALPOSCOEF", type=float, default=0.8)
    parser.add_argument("--torussizecoef", dest="TORUSSIZECOEF", type=float, default=1.2)
    args = parser.parse_args()

    # default values of sizes of petals
    EWIDTH = args.EWIDTH
    EHEIGHT = args.EHEIGHT
    MINWIDTH = args.MINWIDTH
    MINHEIGHT = args.MINHEIGHT

    # shift from the origin(=upper-left end)
    SHIFTX = args.SHIFTX
    SHIFTY = args.SHIFTY

    # coeffcient for the adjustment of petal position
    PETALPOSCOEF = args.PETALPOSCOEF

    # coefficient for the adjustment of torus size
    TORUSSIZECOEF = args.TORUSSIZECOEF

    # coefficients for size calculation
    sizecoef = 0.3

    # read the result file
    significance = -1.0
    try:
        significance = readResult(args.result_filename, args.target_filename, args.value_filename)
    except IOError:
        print(f'Error: "{args.result_filename}" cannot be found.', file=sys.stderr)
        return
    colorvalMin = -math.log(1.0 / significance)

    # for each combination
    for j in range(len(combiName)):

        # reset the end position values
        xmin = 1.0e+30
        xmax = -1.0e+30
        ymin = 1.0e+30
        ymax = -1.0e+30

        # determine the scale of drawing
        nameset = combiName[j]
        for k in range(len(nameset)):
            
            # search for the motif and specify its ID
            nameid = -1
            for i in range(len(motifName)):
                if nameset[k] in motifName[i]:
                    nameid = i
                    break

            # calculate the end position values
            sizeval = math.log(motifNgenes[nameid]) * sizecoef
            rad = k * math.pi * 2 / len(combiName[j]) - math.pi * 0.5
            sizex = EWIDTH * sizeval
            sizey = EHEIGHT * sizeval
            shiftx = (math.cos(rad) * (sizex - EHEIGHT * 2)) + SHIFTX
            shifty = (math.sin(rad) * (sizex - EHEIGHT * 2)) + SHIFTY
            if xmax < shiftx + sizex:
                xmax = shiftx + sizex
            if xmin > shiftx - sizex:
                xmin = shiftx - sizex
            if ymax < shifty + sizey:
                ymax = shifty + sizey
            if ymin > shifty - sizey:
                ymin = shifty - sizey

        # calculate the scaling factor
        scalex = (xmax - xmin) / (SHIFTX * 3)
        scaley = (ymax - ymin) / (SHIFTY * 3)

        # open the SVG file
        nameset = combiName[j]
        svgfilename = f"{args.result_filename}-flower{combiRank[j]}.svg"
        with svg.open_svg(svgfilename) as svgfile:
            # for each motif in the combination
            for k in range(len(nameset)):
                
                # search for the motif and specify its ID
                nameid = -1
                for i in range(len(motifName)):
                    if nameset[k] in motifName[i]:
                        nameid = i
                        break

                # calculate the size and color of the ellipsoid
                sizeval = math.log(motifNgenes[nameid]) * sizecoef
                if motifApvalue[nameid] > significance:
                    colorval = -math.log(motifApvalue[nameid] / significance)
                    if colorval < colorvalMin:
                        colorval = colorvalMin
                if motifApvalue[nameid] <= significance:
                    colorval = motifApvalue[nameid] / significance

                # initialize variables before calculating the shape of the ellipsoid
                rad = k * math.pi * 2 / len(combiName[j]) - math.pi * 0.5
                sizex = EWIDTH * sizeval * scalex
                sizey = EHEIGHT * sizeval * scaley
                if sizex < MINWIDTH:
                    sizex = MINWIDTH
                if sizey < MINHEIGHT:
                    sizey = MINHEIGHT
                shiftx = (math.cos(rad) * (sizex - EHEIGHT * PETALPOSCOEF)) + SHIFTX
                shifty = (math.sin(rad) * (sizex - EHEIGHT * PETALPOSCOEF)) + SHIFTY
                annox = (math.cos(rad) * (sizex * 2  - EHEIGHT * PETALPOSCOEF)) + SHIFTX
                annoy = (math.sin(rad) * (sizex * 2  - EHEIGHT * PETALPOSCOEF)) + SHIFTY

                # draw a motif
                svg.drawMotif(sizex, sizey, shiftx, shifty, rad, colorval, svgfile)
                svg.annotateMotif(motifName[nameid], motifApvalue[nameid],
                                  annox, annoy, svgfile)

            # draw the combination
            colorval = combiApvalue[j] / significance
            size = EHEIGHT * scalex * TORUSSIZECOEF
            if size < MINHEIGHT:
                size = MINHEIGHT
            svg.drawMotif(size, size, SHIFTX, SHIFTY, 0.0, colorval, svgfile)
            svg.annotateMotif("", combiApvalue[j], SHIFTX-30, SHIFTY-5, svgfile)

# call the main function
if __name__ == '__main__':
    main()

