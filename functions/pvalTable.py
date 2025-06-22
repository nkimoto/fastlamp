#!/usr/bin/env python

"""
Copyright (c) 2013, LAMP Team
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * Neither the name of the LAMP Team nor the
      names of its contributors may be used to endorse or promote products
      derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL LAMP TEAM BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""

# Table for storing P-value.
# This source is used in Fisher's exact test and Chi-square test.
# @author Terada, 16, Apr, 2013

import sys

class PvalTable():
	def __init__( self, row_size ):
		self.table = [[-1] * (i + 1) for i in range(row_size)]
	
	def getValue( self, row, col ):
		if len(self.table) <= row:
			return -1
		if len(self.table[row]) <= col:
			return -1
		return self.table[row][col]

	def putValue( self, row, col, pval ):
		if len(self.table) <= row:
			return
		if len(self.table[row]) <= col:
			return
		self.table[row][col] = pval

	def hashSize( self ):
		size = 0
		for row in self.table:
			size = size + len( self.table )
		return size

	def output( self ):
		for i in range( 0, len(self.table) ):
			row = self.table[i]
			sys.stdout.write("[%s]" % i)
			for j in range( 0, len(row) ):
				sys.stdout.write(" %s:%s" % (j, row[j]))
			sys.stdout.write("\n")
	
