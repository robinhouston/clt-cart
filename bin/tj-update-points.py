#!/usr/bin/python
# -*- encoding: utf-8 -*-
from __future__ import division

import json
import optparse
import sys

def update_points(topojson_filename, points_filename, options, fail):
    with open(topojson_filename, 'r') as f:
        data = json.load(f)
    with open(points_filename, 'r') as f:
        points = [
            map(float, line.strip().split(" "))
            for line in f
        ]
    transform = data["transform"]
    translate = transform["translate"]
    scale = transform["scale"]
    arcs = data["arcs"]
    
    points_iter = iter(points)
    new_arcs = []
    for arc in arcs:
        new_arc = []
        new_arcs.append(new_arc)
        px, py = 0, 0
        for (qx, qy) in arc:
            x, y = points_iter.next()
            # Transform to map coordinates
            mx = options.x_min + (x * (options.x_max - options.x_min) / options.width)
            my = options.y_min + (y * (options.y_max - options.y_min) / options.height)
            
            # Transform to topojson grid
            x = (mx - translate[0]) / scale[0]
            y = (my - translate[1]) / scale[1]
            
            new_arc.append([x-px, y-py])
            px, py = x, y
    
    data["arcs"] = new_arcs
    json.dump(data, sys.stdout)

def main():
    parser = optparse.OptionParser(usage="%prog [options] topojson-file points.txt")
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
        parser.error("No topojson or points file specified")
    elif len(args) == 1:
        parser.error("No points file specified")
    elif len(args) > 2:
        parser.error("Unexpected argument: " + args[1])
    
    topojson_filename, points_filename = args
    update_points(topojson_filename, points_filename, options, parser.error)

main()
