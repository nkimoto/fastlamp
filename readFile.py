#!/usr/bin/env python

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

# Define methods to read transaction and flag files
# @author Terada 26, June, 2011
# @editor aika 1, Nov. 2011
#    Developing the readValueFile for handling the value.
# @editor aika 10, Nov. 2011
#    Change readTransactionFile. If file has space the front and the end, remove them.
# @editor aika, 11, Mar. 2014
#    Change readFiles for keeping transaction ID.

import sys
import transaction

##
# Read transaction data and value data.
# This method is simple read.
# Not have any check for input files.
# So, you must prepare complete files.
#
# Return
# 1. list of transactions
# 2. list of mapping from column id to column name
##
def readFiles(transaction_file, value_file, delm):
	gene2id = {} # dictionary from gene name to id
	transaction_list, columnid2name = readItemFile(transaction_file, gene2id, delm)
	transaction_list = readValueFile(value_file, transaction_list, gene2id, delm)
	transaction_list.sort(key=lambda x: x.value)
	return transaction_list, columnid2name

##
# Read item file.
# Return
# 1. list of transactions
# 2. list of mapping from column id to column name
##
def readItemFile(transaction_file, gene2id, delm):
	transaction_list = [] # list of transactions
	columnid2name = [] # list about mapping column id to column name
	gene_set = set([])
	line_num = 0
	col_size = -1
	try:
		with open( transaction_file, 'r', encoding='utf-8' ) as f:
			for line in f:
				line_num = line_num + 1
				row_list = [x.strip() for x in line.strip().split(delm)]
				# If line is header line, read column name
				if line_num == 1:
					col_size = len( row_list )
					for i in range(1, col_size):
						colname = row_list[i]
						columnid2name.append(colname)
					continue
				
				t_name = row_list[0]
				if t_name in gene_set:
					message = "Error: %s is contained two or more times in %s.\\n" % (t_name, transaction_file)
					sys.stderr.write(message)
					raise ValueError(message)
				# check the number of columns
				if len( row_list ) != col_size:
					message = "Error in %s\\n" % transaction_file + \
									 "    The header line contains %s columns, while line %s contains %s columns.\\n" % (col_size, line_num, len( row_list ))
					sys.stderr.write(message)
					raise ValueError(message)
				
				gene_set.add(t_name)
				t = transaction.Transaction(t_name)
				gene2id[t_name] = len(transaction_list)
				for i in range(1, len(row_list)):
					flag = int(row_list[i])
					if flag == 1:
						t.addItem(i)
				transaction_list.append(t)
	except IOError as e:
		sys.stderr.write("Error: %s is not found.\n" % e.filename)
		raise e
	return transaction_list, columnid2name

##
# Read value file.
##
def readValueFile(value_file, transaction_list, gene2id, delm):
	line_num = 0
	gene_set = set([])
	try:
		with open( value_file, 'r', encoding='utf-8' ) as f:
			for line in f:
				line_num = line_num + 1
				if (line.startswith("#")):
					continue
				row_list = [x.strip() for x in line.strip().split(delm)]

				# This error raises if value file contains more than two columns.
				if not len(row_list) == 2:
					message = "Error: line %s in %s.\\n" % (line_num, value_file) + \
									 "       value-file should contain two columns.\\n"
					sys.stderr.write(message)
					raise ValueError(message)

				genename = row_list[0].strip()
				exp_value = row_list[1].strip()

				# This error raises if value cannot be converted to float.
				try:
					exp_value = float(exp_value)
				except ValueError as e:
					message = "Error: line %s in %s.\\n" % (line_num, value_file) + \
									 "       %s cannot be converted to float.\\n" % (exp_value)
					sys.stderr.write(message)
					raise e

				gene_set.add(genename)
				transaction_list[ gene2id[genename] ].setValue(exp_value)
	except IOError as e:
		sys.stderr.write("Error: %s cannot be opened.\n" % e.filename)
		raise e
	return transaction_list

##
# Check two files transaction name are same
# If transaction contains only one file, flag is -1.
##
def checkTransName(transaction_list, transaction_file):
	for t in transaction_list:
		if t.value is None:
			message = "\"%s\" only appears in %s\n" % (t.name, transaction_file)
			sys.stderr.write(message)
			raise ValueError(message)


##
#
##
def colname2id(columnid2name):
	colname2id_dict = {}
	index = 0
	for i in columnid2name:
		colname2id_dict[i] = index
		index = index + 1
	return colname2id_dict

def isFloat( value_str ):
	try:
		float( value_str )
		return True
	except ValueError:
		return False
