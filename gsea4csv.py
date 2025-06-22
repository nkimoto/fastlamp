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

# Convert the C3 of GSEA to csv format.
# The made file represents the association between gene and TF.
# "," is deliminated in the output file.
# The first column is a gene ID of Entrez gene.
# From the second columns,  the value indicates whether the gene is regulated the TF.
# If the gene is regulated the TF, then the value is 1, otherwise 0.
# @author aika, 16, Nov., 2011

__author__ = "Aika TERADA"

import sys
import argparse


##
# read GMT file and return gene dictionary.
# The dictionary map gene name and interacted TFs list.
##
def readGmtFile( gmt_file ):
    gene_dict = {} # key: gene, value: TFs list
    all_motif_lst = [] # TFs list which appear in gme file
    with open(gmt_file, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            s = line.split("\t")
            tf_set = s[0] # TFs set include represent motif
            motif = ""
            # If the TF does not known
            if "$" not in tf_set:
                motif = tf_set
            else:
                s_tf_set = tf_set.split("$")
                motif = s_tf_set[1].split("_")[0]

            if motif not in all_motif_lst:
                all_motif_lst.append(motif)

            # mapping gene and the tf motif
            interacted_genes = s[2:] # interacted genes
            for gene in interacted_genes:
                if gene not in gene_dict:
                    gene_dict[gene] = []
                if motif not in gene_dict[gene]:
                    gene_dict[gene].append(motif)
    return all_motif_lst, gene_dict

##
# output CSV format.
# all_tf_list: List of all TFs included GMT file.
# gene_dict: [key] gene, [value] regulate TF list.
# output_file: The filename to print CSV format.
##
def outCSVFormat( all_tf_lst, gene_dict, output_file ):
    try:
        with open( output_file, "w", encoding='utf-8' ) as f:
            # output header
            f.write( "#EntrezGene" )
            for tf in all_tf_lst:
                f.write(f",{tf}")
            f.write("\n")

            for gene in gene_dict:
                f.write(f"{gene}")
                regulated_tfs = gene_dict[ gene ]
                for tf in all_tf_lst:
                    if tf in regulated_tfs:
                        f.write( ",1" )
                    else:
                        f.write( ",0" )
                f.write("\n")
    except OSError as e:
        sys.stderr.write("Error in output CSV file.\n")
        raise e

def run( gmt_file, output_file ):
    print( "Read GMT file...", file=sys.stderr )
    all_tf_list, gene_dict = readGmtFile( gmt_file )
    print( "Make GSV file...", file=sys.stderr )
    outCSVFormat(all_tf_list, gene_dict, output_file)
    print(f"# of TFs: {len(all_tf_list)}, # of genes: {len(gene_dict)}")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description="Convert the C3 of GSEA to csv format.")
    p.add_argument("gmt_file", help="input gmt file")
    p.add_argument("output_file", help="output csv file")

    args = p.parse_args()

    try:
        run( args.gmt_file, args.output_file )
    except (OSError, IOError) as e:
        sys.stderr.write(f"Error: {e}\n")
        sys.exit(1)
