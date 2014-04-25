#!/usr/bin/python
# -*- encoding: utf-8 -*-
from __future__ import division

import csv
import math
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
    
    sum_of_values = 0
    sum_of_squares = 0
    n = 0
    
    for k in input_data.iterkeys():
        if k in area_data:
            area = area_data[k]
            density = area / input_data[k]
            sum_of_values += density
            sum_of_squares += density*density
            n += 1
    
    mean = sum_of_values / n
    variance = (sum_of_squares / n) - mean*mean
    print "Mean = %g, standard deviation = %g, n=%d" % (mean, math.sqrt(variance), n)

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
