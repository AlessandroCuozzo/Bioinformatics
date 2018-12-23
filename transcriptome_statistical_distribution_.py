#!/usr/bin/env python3
#-*- coding: utf-8 -*-

"""
Just the basics of a transcriptome FASTA file statistical distribution :
- number of contigs
- minimum, maximum, average and median lengths
- N50 and L50

copyright (c) 2018 Alessandro Cuozzo

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

"""

# https://en.wikipedia.org/wiki/N50,_L50,_and_related_statistics
# http://www.acgt.me/blog/2015/6/11/l50-vs-n50-thats-another-fine-mess-that-bioinformatics-got-us-into

# ================================
# READ FILE
# ================================

from numpy import*

filename="name_of_your_fasta_file"
dico={}

with open(filename,"r") as f:
	for line in f:
		new_line=line.rstrip()
		if new_line[0]==">":
			header=new_line
			dico[header]=""
		else:
			dico[header]=dico.get(header,"")+new_line


# ================================
# Basic distribution
# ================================

liste_len=[]
for N in dico.values():
	liste_len.append(float(len(N)))

contigs = len(liste_len)
average = sum(liste_len)/contigs
MIN = min(liste_len)
MAX = max(liste_len)
median = sorted(liste_len)[(contigs/2)-1]

# ================================
# N50 & L50
# ================================

T = sum(liste_len)/2 # Threshold
short = sorted(list(set(liste_len))) #[::-1]
M = array(sorted(liste_len) #[::-1])

N50 = 0 # LENGTH of contigs >= 50% genome size
for L in short: # L = Length
	if M[where(M>=L)].sum() >= T and L > N50:
		N50 = L

L50 = len(M[M-N50>=0]) # NUMBER of contigs >= 50% genome size

# =====
# print
# =====

print("Ncontigs\taverage\tmin\tmax\tmedian\tN50\tL50")
print(contigs,'\t',average,'\t',MIN,'\t',MAX,'\t',median,'\t',N50,'\t',L50)
