#!/usr/bin/env python3
#-*- coding: utf-8 -*-

# ===================================================================
# Translation of RNA sequences into protein sequences ||| STOP === X
# ===================================================================

"""This translation aims to keep track of every potential peptide encoded in a transcript. Insetad of selecting only the longest in silico protein, we will write all peptides found, for all six frameshifts 3*(5'->3') + 3*(3'->5' reverse complement).

The purpose of this methodology is to consider an overlapping translation ; we assume that a single transcript can encode multiple peptides.

It uses the biopython module to perform the translation.

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

code={'TTT': 'F','TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'TCT': 'S', 'TCC': 'S', 'TCA': 'S','TCG': 'S','TAT': 'Y','TAC': 'Y','TAA': 'X','TAG': 'X','TGT': 'C','TGC': 'C','TGA': 'X','TGG': 'W','CTT': 'L','CTC': 'L','CTA': 'L','CTG': 'L','CCT': 'P','CCC': 'P','CCA': 'P','CCG': 'P','CAT': 'H','CAC': 'H','CAA': 'Q','CAG': 'Q','CGT': 'R','CGC': 'R','CGA': 'R','CGG': 'R','ATT': 'I','ATC': 'I','ATA': 'I','ATG': 'M','ACT': 'T','ACC': 'T','ACA': 'T','ACG': 'T','AAT': 'N','AAC': 'N','AAA': 'K','AAG': 'K','AGT': 'S','AGC': 'S','AGA': 'R','AGG': 'R','GTT': 'V','GTC': 'V','GTA': 'V','GTG': 'V','GCT': 'A','GCC': 'A','GCA': 'A','GCG': 'A','GAT': 'D','GAC': 'D','GAA': 'E','GAG': 'E','GGT': 'G','GGC': 'G','GGA': 'G','GGG': 'G'}

from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

filename = "FASTA_file_to_read"
filename2 = "new_FASTA_file_to_write"

# ========================================
# BUILD DICTIONNARY
# ========================================

dico={}

with open(filename,"r") as f:
	for line in f:
		new_line=line.rstrip()
		if len(new_line)>0:
			if new_line[0]==">":
				header=new_line
				dico[header]=""
			else:
				dico[header]=dico.get(header,"")+new_line


# ========================================
# TRANSLATION + WRITE IN FILE
# ========================================


liste=[]

with open(filename2, "w") as f2:
	for sequence in sorted(dico.keys()):
		proteine=""
		seq=dico[sequence].upper().replace('U','T')  # 5' -> 3'
		frames=[]
		for n in range(2):  # Once for 5' -> 3' and once for 3' -> 5' reverse complement
			for frameshift in range(3):
				while len(seq)%3 !=0:
					seq=seq[0:-1]
				coding_dna = Seq(seq,IUPAC.unambiguous_dna)
				pep = str(coding_dna.translate())
				seq=seq[1:]
				frames.append(pep)

			if n==1:
				break
			seq=dico[sequence][::-1].upper().replace('U','T') # 3'-> 5' / reverse complement
			seq2=""
			for element in seq:
				if element=="A":
					seq2+="T"
				elif element=="T":
					seq2+="A"
				elif element=="C":
					seq2+="G"
				elif element=="G":
					seq2+="C"
			seq=seq2

		E=0 # frame (there will be six for each query contig [transcript])
		for proteine in frames:
			E+=1
			f2.write(sequence.split("_")[0]+"_"+str(E)+"|"+sequence.split("_")[1]+"\n"+proteine+"\n")




