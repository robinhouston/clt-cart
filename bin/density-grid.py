#!/usr/bin/python
# -*- encoding: utf-8 -*-
from __future__ import division

import csv
import math
import optparse

NaN = float("NaN")

def load_csv(filename, key_col, value_col, fail, value_type):
    ret = {}
    with open(filename, 'r') as f:
        r = csv.reader(f)
        header = r.next()
        if key_col not in header:
            fail("File '%s' does not have a column '%s'" % (filename, key_col))
        if value_col not in header:
            fail("File '%s' does not have a column '%s'" % (filename, value_col))
        for row in r:
            d = dict(zip(header, row))
            ret[d[key_col]] = value_type(d[value_col])
    return ret

def density_grid(grid_filename, numerator_filename, denominator_filename, options, fail):
    with open(grid_filename, 'r') as f:
        region_grid = list(csv.reader(f))
    numerators = load_csv(numerator_filename, options.numerator_key, options.numerator_value, fail, float)
    denominators = load_csv(denominator_filename, options.denominator_key, options.denominator_value, fail, float)
    
    grid = [
        [
            numerators.get(region, NaN) / denominators.get(region, 1.0)
            for region in row[1:]
        ]
        for row in region_grid
    ]
    
    n, total = 0, 0.0
    for row in grid:
        for value in row:
            if not math.isnan(value):
                n += 1
                total += value
    mean_value = total / n
    
    return [
        [
            1.0 if math.isnan(value) else value / mean_value
            for value in row
        ]
        for row in grid
    ]

def main():
    parser = optparse.OptionParser(usage="%prog [options] grid.csv numerator.csv denominator.csv")
    parser.add_option("", "--numerator-key",
                    action="store", type="str",
                    default="key",
                    help="Column name of numerator key (default: %default)")
    parser.add_option("", "--numerator-value",
                    action="store", type="str",
                    default="value",
                    help="Column name of numerator value (default: %default)")
    parser.add_option("", "--denominator-key",
                    action="store", type="str",
                    default="key",
                    help="Column name of denominator key (default: %default)")
    parser.add_option("", "--denominator-value",
                    action="store", type="str",
                    default="value",
                    help="Column name of denominator value (default: %default)")
    (options, args) = parser.parse_args()
    
    if len(args) < 3:
        parser.error("Not enough arguments")
    elif len(args) > 3:
        parser.error("Unexpected argument: " + args[3])
    
    grid_filename, numerator_filename, denominator_filename = args
    grid = density_grid(grid_filename, numerator_filename, denominator_filename,
        options, parser.error)
    
    width, height = len(grid[0]), len(grid)
    for row in grid:
        print " ".join(map(str, row))

main()
