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

# Define fuctions that is used in multiple_test.py
# Calculate the P-value and lower bound of the combination. 
# @author Terada, 26, June, 2011
# @editor Terada, 24, Feb, 2012,
#     If x > n1 in caluclation MASL, then raise error and exit (This case does not treat the case).
# @editor Terada, 16, Apr, 2013,
#     Acceletate of the calculation P-value by storing the calculated P-value.
# @editor Terada, 11, Mar, 2015,
#     Implement computation of the 'less' and 'two-sided' Fisher's exact test. 
from __future__ import division
import sys
import os
from . import functionsSuper as fs
from . import pvalTable
from scipy.stats import hypergeom
from numba import njit
import math

pardir = os.path.dirname(os.path.dirname(os.path.abspath( __file__ )))
sys.path.append(pardir)

# Numba-optimized functions for Fisher's exact test
@njit(cache=True)
def _numba_log_factorial(n):
    """Numba-optimized log factorial calculation"""
    if n <= 1:
        return 0.0
    result = 0.0
    for i in range(2, n + 1):
        result += math.log(i)
    return result

@njit(cache=True)
def _numba_log_combination(n, k):
    """Numba-optimized log combination calculation: log(C(n,k))"""
    if k > n or k < 0:
        return float('-inf')
    if k == 0 or k == n:
        return 0.0
    if k > n - k:
        k = n - k  # Use symmetry for efficiency
    
    result = 0.0
    for i in range(k):
        result += math.log(n - i) - math.log(i + 1)
    return result

@njit(cache=True)
def _numba_hypergeom_pmf(k, n, K, N):
    """Numba-optimized hypergeometric PMF calculation"""
    # P(X = k) = C(K, k) * C(N-K, n-k) / C(N, n)
    if k < max(0, n + K - N) or k > min(n, K):
        return 0.0
    
    log_prob = (_numba_log_combination(K, k) + 
               _numba_log_combination(N - K, n - k) - 
               _numba_log_combination(N, n))
    
    return math.exp(log_prob)

@njit(cache=True)
def _numba_hypergeom_sf(k, n, K, N):
    """Numba-optimized hypergeometric survival function: P(X >= k)"""
    if k < 0:
        return 1.0
    
    max_k = min(n, K)
    min_k = max(0, n + K - N)
    
    if k > max_k:
        return 0.0
    if k <= min_k:
        return 1.0
    
    total_prob = 0.0
    for i in range(k, max_k + 1):
        total_prob += _numba_hypergeom_pmf(i, n, K, N)
    
    return total_prob

@njit(cache=True)
def _numba_fisher_exact_greater(a, b, c, d):
    """Numba-optimized Fisher's exact test (greater)"""
    n = a + b  # row total
    K = a + c  # column total  
    N = a + b + c + d  # grand total
    return _numba_hypergeom_sf(a, n, K, N)

@njit(cache=True)
def _numba_fisher_exact_two_sided(a, b, c, d):
    """Numba-optimized Fisher's exact test (two-sided)"""
    n = a + b
    K = a + c
    N = a + b + c + d
    
    # Get probability of observed value
    observed_prob = _numba_hypergeom_pmf(a, n, K, N)
    
    # Sum all probabilities <= observed probability
    total_prob = 0.0
    max_k = min(n, K)
    min_k = max(0, n + K - N)
    
    for k in range(min_k, max_k + 1):
        prob_k = _numba_hypergeom_pmf(k, n, K, N)
        if prob_k <= observed_prob + 1e-12:  # Small tolerance for floating point
            total_prob += prob_k
    
    return min(total_prob, 1.0)


##
# Define class
# This class calculate function f that define in paper. (C(n1, x)/C(n0+n1, x))
# transaction_list: list of transactions
##
class FunctionOfX(fs.FunctionsSuper):
    def __init__(self, transaction_list, row_size, alternative):
        fs.FunctionsSuper.__init__(self)
        self.__t_size = len(transaction_list) # all transaction size
        self.__f_size = self.sumValue(transaction_list) # transaction size which have flag = 1 (n1)
        self.__pvalTable = pvalTable.PvalTable( row_size ) # P-value table
        self.__occrTable = pvalTable.PvalTable( row_size ) # occurence table for calculate P-value
        self.calTime = 0 # Total number of calculate P-value
        self.alternative = alternative # alternative hypothesis. greater or less -> 1, two.sided -> 0.
        if self.__f_size == 0:
            message = "Error: There is no up-regulate gene.\n"
            sys.stdout.write(message)
            raise ValueError(message)
        # Check the transaction value.
        # If the value is not 1 or 0, raise error.
        # Because fisher's exact test does not handle numerical value.
        for t in transaction_list:
            if not (t.value == 1.0 or t.value == 0.0):
                message = "Error: \"" + t.name + "\" value is " + str(t.value)+".\n" + \
									"       But value is 1 or 0 if you test by fisher's exact test.\n"
                sys.stderr.write(message)
                raise ValueError(message)

            
    def getN1(self):
        return self.__f_size
    
    def getAllSize(self):
        return self.__t_size

    
    def funcF(self, x):
        """High-precision MASL calculation for Fisher exact test"""
        M = self.__t_size
        n = self.__f_size
        N = x
        k = x
        
        # Ensure we use scipy for highest precision in MASL calculations
        return hypergeom.sf(k - 1, M, n, N)
    
            
    ##
    # Calculate p-value by using fisher's exact test.
    # transaction_list: List of transactions
    # flag_itemset_id: Transactions which have items
    ##
    def calPValue(self, transaction_list, flag_transactions_id):
        ovalues = self.contingencyTable(transaction_list, flag_transactions_id, self.__t_size, self.__f_size)
        a, b, c, d = ovalues[0][0], ovalues[0][1], ovalues[1][0], ovalues[1][1]
        
        # Always use Numba-optimized version
        if self.alternative == 0:  # two-sided
            p_value = _numba_fisher_exact_two_sided(a, b, c, d)
        else:  # greater
            p_value = _numba_fisher_exact_greater(a, b, c, d)
        
        return p_value, a
    
    def _ultra_fast_hypergeom_sf(self, x, n, n1, k):
        """Ultra-fast survival function for hypergeometric distribution: P(X >= x)"""
        import math
        
        # Ultra-fast edge cases
        if x <= 0:
            return 1.0
        
        max_possible = min(n1, k)
        min_possible = max(0, n1 + k - n)
        
        if x < min_possible:
            return 1.0
        if x > max_possible:
            return 0.0
        if x == max_possible:
            # Only one possible value, calculate directly
            try:
                log_prob = (math.lgamma(n1 + 1) - math.lgamma(x + 1) - math.lgamma(n1 - x + 1) +
                           math.lgamma(n - n1 + 1) - math.lgamma(k - x + 1) - math.lgamma(n - n1 - k + x + 1) -
                           (math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)))
                return min(1.0, math.exp(log_prob))
            except (OverflowError, ValueError):
                return hypergeom.sf(x - 1, n, n1, k)
        
        # Use efficient calculation for small ranges
        range_size = max_possible - x + 1
        if range_size <= 3:
            # Direct calculation for very small ranges
            total_prob = 0.0
            for i in range(x, max_possible + 1):
                try:
                    log_prob = (math.lgamma(n1 + 1) - math.lgamma(i + 1) - math.lgamma(n1 - i + 1) +
                               math.lgamma(n - n1 + 1) - math.lgamma(k - i + 1) - math.lgamma(n - n1 - k + i + 1) -
                               (math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)))
                    total_prob += math.exp(log_prob)
                except (OverflowError, ValueError):
                    return hypergeom.sf(x - 1, n, n1, k)
            return min(1.0, total_prob)
        
        # Fallback to scipy for complex cases
        return hypergeom.sf(x - 1, n, n1, k)
    
    def _ultra_fast_hypergeom_two_sided(self, x, n, n1, k):
        """Ultra-fast two-sided test for hypergeometric distribution"""
        try:
            import math
            
            # Ultra-fast edge cases
            max_possible = min(n1, k)
            min_possible = max(0, n1 + k - n)
            
            if max_possible == min_possible:
                return 1.0  # Only one possible value
            
            # Get the probability of the observed value
            log_prob_observed = (math.lgamma(n1 + 1) - math.lgamma(x + 1) - math.lgamma(n1 - x + 1) +
                               math.lgamma(n - n1 + 1) - math.lgamma(k - x + 1) - math.lgamma(n - n1 - k + x + 1) -
                               (math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)))
            prob_observed = math.exp(log_prob_observed)
            
            # For small ranges, use direct calculation
            range_size = max_possible - min_possible + 1
            if range_size <= 10:
                total_prob = 0.0
                for i in range(min_possible, max_possible + 1):
                    log_prob_i = (math.lgamma(n1 + 1) - math.lgamma(i + 1) - math.lgamma(n1 - i + 1) +
                                 math.lgamma(n - n1 + 1) - math.lgamma(k - i + 1) - math.lgamma(n - n1 - k + i + 1) -
                                 (math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)))
                    prob_i = math.exp(log_prob_i)
                    
                    # Include if probability is less than or equal to observed (more extreme)
                    if prob_i <= prob_observed + 1e-12:  # Small tolerance for floating point comparison
                        total_prob += prob_i
                
                return min(total_prob, 1.0)
            
            # Fallback to scipy for complex cases
            return hypergeom.sf(x - 1, n, n1, k)
            
        except Exception:
            # Final fallback to scipy
            return hypergeom.sf(x - 1, n, n1, k)
    
    ##
    # Calculate probability of occurrence probability about table.
    # x: total of the first row (the number of targetting gene)
    # a: the top-left of the table
    # Return C(n1, a)*C(n0, b)/C(n1+n0, a+b)
    ##
    def __probability(self, x, a):
        p = self.__occrTable.getValue(x, a)
        if p < 0:
            n = self.__t_size
            n1 = self.__f_size
            p = hypergeom.pmf(a, n, n1, x)
            self.__occrTable.putValue(x, a, p)
        return p

def run(xls_file, value_file, itemset_str_lst, delimiter, alternative):
    global readFile
    import readFile
    transaction_list, columnid2name = readFile.readFiles(xls_file, value_file, delimiter)
    max_lambda = maxLambda(transaction_list)
    if alternative < 0:
        global lamp
        from lamp import reverseValue
        transaction_list = reverseValue( transaction_list, "fisher" )
    func = FunctionOfX(transaction_list, max_lambda, abs( alternative ) )
    colname2id_dict = readFile.colname2id(columnid2name)

    itemset = set()
    for i in itemset_str_lst:
        item_id = colname2id_dict[i]
        itemset.add(item_id + 1)
        
    flag_transactions_id = []
    for i in range( len(transaction_list) ):
        t = transaction_list[i]
        if len( itemset & t.itemset ) == len(itemset):
            flag_transactions_id.append( i )
    p_value, stat_value = func.calPValue(transaction_list, flag_transactions_id)
    n = len(transaction_list)
    n1 = func.getN1()
    sys.stdout.write("p-value: %s (N: %s, n1: %s, x: %s, a: %s)\n"
                     % (p_value, n, n1, len(flag_transactions_id), stat_value))
    return (p_value, len(flag_transactions_id))

##
# Return max lambda. That is, max size itemset.
##
def maxLambda(transaction_list):
    # Count each item size
    item_sizes = {}
    for t in transaction_list:
        for item in t.itemset:
            # If item does not exist in item_size, then make mapping to 0
            if item not in item_sizes:
                item_sizes[item] = 0
            item_sizes[item] = item_sizes[item] + 1
    
    # Get max value in item_sizes
    max_value = 0
        # itervalues of dict was depricated
        # It is replaced as values for getting each of value in item_sizes
    # for i in item_sizes.itervalues():
    for i in item_sizes.values():
        if i > max_value:
            max_value = i
            
    return max_value


if __name__ == "__main__":
    if (len(sys.argv) < 4):
        sys.stderr.write("Error: functions4fisher.py [item-file] [value-file] [itemset]\n")
        sys.exit()

    xls_file = sys.argv[1]
    value_file = sys.argv[2]
    itemset_str_lst = sys.argv[3].split(',')
    delimiter = ','
    p_value, down_size = run(xls_file, value_file, itemset_str_lst, delimiter)
