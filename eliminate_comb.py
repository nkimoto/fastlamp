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

# Remove the redundant combinations in the result of LAMP.
# @author Terada 3, Aug., 2013


import sys
import argparse
from operator import itemgetter

# read the result of LAMP.
# filename: the result filename of LAMP
def readResult( filename ):
    try:
        detections_list = [] # Store the combinations 
        meta_line_list = [] # Store the meta lines
        time_line = "" # Store the line containing running time
        minp_line_list = [] # Store the minimum P-value lines
        with open(filename, 'r', encoding='utf-8') as f:
            flag_broken = True # a flag to decide whether the file is broken or not.
            flag_comb = False
            flag_minp_dist = False # a flag for checking that the line is minimum P-value distribution.
            for line in f:
                line = line.strip()
                if ( line.startswith("Time (sec.): ") ):
                    time_line = line
                    flag_broken = False
                    flag_comb = False
                    continue
                if line == "--- minimum P-values ---":
                    minp_line_list.append( line )
                    flag_minp_dist = True
                    continue
                if flag_minp_dist:
                    minp_line_list.append( line )
                    continue

                if not (flag_comb or flag_minp_dist):
                    meta_line_list.append( line )
                    if ( line.startswith("Rank") ):
                        flag_comb = True
                    continue
                # [0]: Rank, [1]: Raw p-value, [2]: Adjusted P-value,
                # [3]: detections, [4]: arity, [5]: support, [6]: statistic_score.
                detections = line.split('\t')
                detections[1] = float(detections[1])
#                detections[2] = float(detections[2])
                detections[4] = int(detections[4])
                detections_list.append( tuple( detections ) )
            # if the file does not have results, output error.
            if (flag_broken):
                sys.stderr.write("%s is broken.\n" % filename)
        return detections_list, meta_line_list, time_line, minp_line_list
    except IOError:
        message = "'%s' cannot be opened.\n" % filename
        sys.stderr.write(message)
        raise IOError(message)

# Return whether set_j is subset of set_i.
def isSubset(set_i, set_j):
#    print set_i
#    print set_j
    diff_set = set_i - set_j
    if (len(diff_set) == (len(set_i) - len(set_j))):
#        print "subset"
        return True
    else:
        return False

# Eliminate the redundant combinations.
# detections_list: a list of combinations made by readResult function.
def mergeResult( detections_list ):
    upper_list = []
    merged_list = []
    for detection in detections_list:
        detect_set = set(detection[3].split(','))
        # check whether the detection is subset of the list which is already found.
        flag = True
        for i in upper_list:
            if ( isSubset(detect_set, i) ) or ( isSubset(i, detect_set) ):
                flag = False
                break
        if flag:
            merged_list.append( detection )
        upper_list.append( detect_set )
    return merged_list

# Sort detections_list by P-value.
# If comb1 and comb2 have equal P-values, the larger combination gives high rank.
def sortComb( detections_list ):
    detections_list = sorted( detections_list, key = itemgetter( 1, 4 ) )
    return detections_list

# output result
def output( out_filename, detections_list, meta_line_list, time_line, minp_line_list ):
    output_file = sys.stdout
    if len( out_filename ) > 0:
        try:
            output_file = open(out_filename, 'w', encoding='utf-8')
        except IOError:
            message = "Error: Cannot output to %s'\n" % out_filename
            sys.stderr.write(message)
            raise IOError(message)

    print("# Non-redundant combinations", file=output_file)
    for i in meta_line_list:
        print(i, end="", file=output_file)
        if i.startswith("# # of significant combinations: "):
            print(" -> %d" % len(detections_list), end="", file=output_file)
        print(file=output_file)
    rank = 0
    for detect in detections_list:
        rank = rank + 1
        print(rank, end="", file=output_file)
        for i in detect[1:]:
            print("\t%s" % i, end="", file=output_file)
        print(file=output_file)
    print(time_line, file=output_file)
    for i in minp_line_list:
        print(i, file=output_file)

    if output_file is not sys.stdout:
        output_file.close()

# in_filename: the filename of LAMO
# out_filename: the filename to output the eliminated result.
#               When out_filename is "", the result is printed as standard output.
def run( in_filename, out_filename ):
    detections_list, meta_line_list, time_line, minp_line_list = readResult( in_filename )
    detections_list = sortComb( detections_list )
    merged_list = mergeResult( detections_list )
    output( out_filename, merged_list, meta_line_list, time_line, minp_line_list )
    

if __name__ == "__main__":
    p = argparse.ArgumentParser(description="Remove the redundant combinations in the result of LAMP.")
    p.add_argument("filename", help="The result filename of LAMP")
    p.add_argument('-o', '--output', dest = "output_filename", default = "", help = "The file name to output the eliminated result.")
    
    args = p.parse_args()
    
    try:
        run( args.filename, args.output_filename )
    except (ValueError, IOError) as e:
        sys.stderr.write("Error: %s\n" % e)
        sys.exit(1)
