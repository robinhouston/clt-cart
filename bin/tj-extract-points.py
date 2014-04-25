#!/usr/bin/python
# -*- encoding: utf-8 -*-
from __future__ import division

import json
import optparse
import sys

def extract_points(filename, options, fail):
    with open(filename, 'r') as f:
        data = json.load(f)
    transform = data["transform"]
    translate = transform["translate"]
    scale = transform["scale"]
    arcs = data["arcs"]
    
    for arc in arcs:
        px, py = 0, 0
        for (qx, qy) in arc:
            x = (px + qx) * scale[0] + translate[0]
            y = (py + qy) * scale[1] + translate[1]
            
            # Transform to grid coordinates
            gx = (x - options.x_min) * options.width / (options.x_max - options.x_min)
            gy = (y - options.y_min) * options.height / (options.y_max - options.y_min)
            print gx, gy
            
            px, py = px + qx, py + qy

def main():
    parser = optparse.OptionParser(usage="%prog [options] topojson-file")
    parser.add_option("", "--width",
                    action="store", type="int",
                    default=500,
                    help="grid width (default %default)")
    parser.add_option("", "--height",
                    action="store", type="int",
                    default=250,
                    help="grid height (default %default)")
    
    parser.add_option("", "--x-min",
                    action="store", type="float",
                    default=-17005833.33053,
                    help="min x value in input (default %default)")
    parser.add_option("", "--x-max",
                    action="store", type="float",
                    default=17005833.33053,
                    help="max x value in input (default %default)")
    parser.add_option("", "--y-min",
                    action="store", type="float",
                    default=-8625154.47185,
                    help="min x value in input (default %default)")
    parser.add_option("", "--y-max",
                    action="store", type="float",
                    default=8625154.47185,
                    help="max x value in input (default %default)")
    
    (options, args) = parser.parse_args()
    
    if not args:
        parser.error("No topojson file specified")
    elif len(args) != 1:
        parser.error("Unexpected argument: " + args[1])
    
    filename ,= args
    extract_points(filename, options, parser.error)

main()
