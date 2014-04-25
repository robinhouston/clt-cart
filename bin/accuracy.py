#!/usr/bin/python
# -*- encoding: utf-8 -*-
from __future__ import division

import csv
import optparse
import sys

def load(filename, key_name, value_name, fail):
    with open(filename, 'r') as f:
        r = csv.reader(f)
        header = r.next()
        if key_name not in header:
            fail("File %s does not have a column labelled '%s'", filename, key_name)
        if value_name not in header:
            fail("File %s does not have a column labelled '%s'", filename, value_name)
        ret = {}
        for row in r:
            d = dict(zip(header, row))
            ret[d[key_name]] = float(d[value_name])
        
        return ret

def compute_errors(input_filename, areas_filename, options, fail):
    input_data = load(input_filename, options.input_key, options.input_value, fail)
    area_data = load(areas_filename, options.area_key, options.area_value, fail)
    
    min_density = float("Inf")
    max_density = 0
    min_region, max_region = None, None
    
    for k in input_data.iterkeys():
        if k in area_data:
            density = input_data[k] / area_data.get(k)
            if density < min_density:
                min_region = k
                min_density = density
            if density > max_density:
                max_region = k
                max_density = density
    
    print "Density ranges from %g (in %s) to %g (in %s)" % (min_density, min_region, max_density, max_region)

def main():
    parser = optparse.OptionParser(usage="%prog [options] input-data.csv areas.csv")
    parser.add_option("", "--input-key",
                    action="store", type="str",
                    default="Country",
                    help="Key column in input file (default %default)")
    parser.add_option("", "--input-value",
                    action="store", type="str",
                    default="Population",
                    help="Value column in input file (default %default)")
    
    parser.add_option("", "--area-key",
                    action="store", type="str",
                    default="Feature ID",
                    help="Key column in area file (default %default)")
    parser.add_option("", "--area-value",
                    action="store", type="str",
                    default="Area",
                    help="Value column in area file (default %default)")
    
    (options, args) = parser.parse_args()
    
    if not args:
        parser.error("No files specified")
    elif len(args) == 1:
        parser.error("Only one input file specified")
    elif len(args) > 2:
        parser.error("Unexpected argument: " + args[2])
    
    input_filename, areas_filename ,= args
    compute_errors(input_filename, areas_filename, options, parser.error)

main()
