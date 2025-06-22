#!/usr/bin/env python3

"""Copyright (c) 2013, LAMP development team
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

# Compare gene-ID between expression-file and association-file
# and make new expression-file and association-file which are consisted from same gene ID.
# The gene IDs are based on the input expression-file.

import os
import sys
import argparse

__author__ = "Aika TERADA"

SEPARATOR_EXP = ","
SEPARATOR_CSV = ","

##
# read expression-file and return the set of gene IDs in expression-file
# exp_file: The expression filename ( inputted filename )
##
def readEXPFile( exp_file ):
    gene_list = []
    try:
        with open( exp_file, 'r', encoding='utf-8' ) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                s = line.strip().split(SEPARATOR_EXP)
                gene_list.append(s[0] )
        return gene_list
    except OSError as e:
        sys.stderr.write( "Error in read expression-file.\n" )
        raise e

def readCSVFile( csv_file ):
    association_dict = {}
    try:
        with open( csv_file, 'r', encoding='utf-8' ) as f:
            column_line = ""
            for line in f:
                line = line.strip()
                if line.startswith("#"):
                    column_line = line
                    continue
                gene = line.split( SEPARATOR_CSV )[0]
                association_dict[ gene ] = line
        return column_line, association_dict
    except OSError as e:
        sys.stderr.write( "Error in read csv-file.\n" )
        raise e

##
# read CSV file and output a new file when the gene included gene_set
##
def makeCSVFile( out_csv_file, gene_list, association_dict, column_line ):
    tf_size = len( column_line.split(SEPARATOR_CSV) ) - 1
    try:
        with open( out_csv_file, "w", encoding='utf-8' ) as fo:
            fo.write(f"{column_line}\n")
            for gene in gene_list:
                if gene in association_dict:
                    fo.write(f"{association_dict[gene]}\n")
                else:
                    fo.write(f"{gene}")
                    fo.write(",0" * tf_size)
                    fo.write("\n")
    except OSError as e:
        sys.stderr.write( "Error in make a new csv file.\n" )
        raise e

def run( exp_file, csv_file,  out_csv_file):
    gene_list = readEXPFile( exp_file )
    column_line, association_dict = readCSVFile( csv_file )
    makeCSVFile( out_csv_file, gene_list, association_dict, column_line )
    print(f"# of TFs: {len(column_line.split(','))-1}, # of genes : {len(association_dict) - 1}")

if __name__ == "__main__":
    p = argparse.ArgumentParser(description="Compare gene-ID between expression-file and association-file and make new expression-file and association-file which are consisted from same gene ID.")
    p.add_argument("expression_file")
    p.add_argument("csv_file")
    p.add_argument("output_csv_file")
    
    args = p.parse_args()

    # Check the file existence.
    if not os.path.isfile(args.expression_file):
        sys.stderr.write(f"IOError: No such file: '{args.expression_file}'\n")
        sys.exit(1)
    if not os.path.isfile(args.csv_file):
        sys.stderr.write(f"IOError: No such file: '{args.csv_file}'\n")
        sys.exit(1)

    try:
        run( args.expression_file, args.csv_file, args.output_csv_file )
    except (OSError, IOError) as e:
        sys.stderr.write(f"Error: {e}\n")
        sys.exit(1)
