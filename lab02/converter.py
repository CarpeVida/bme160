#!/usr/bin/env python3
# Name: Daniel Schmelter (dschmelt)
# Group Members: None

"""
This program should take in input of a codon or amino acid 1 letter or 3 letter code
and then output the corresponding codon or amino acid matching it via dictionaries
"""
import sys

class Converter:

    #long amino acid code keys, short values
    short_AA = {
                'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
                'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
                'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
                'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'
                }
    #switches keys and values so short key has long AA value
    long_AA = {value:key for key,value in short_AA.items()}

    rnaCodonTable = {
    # Second Base
    # U             C             A             G
    #U
    'UUU': 'Phe', 'UCU': 'Ser', 'UAU': 'Tyr', 'UGU': 'Cys',
    'UUC': 'Phe', 'UCC': 'Ser', 'UAC': 'Tyr', 'UGC': 'Cys',
    'UUA': 'Leu', 'UCA': 'Ser', 'UAA': '---', 'UGA': '---',
    'UUG': 'Leu', 'UCG': 'Ser', 'UAG': '---', 'UGG': 'Trp',
    #C
    'CUU': 'Leu', 'CCU': 'Pro', 'CAU': 'His', 'CGU': 'Arg',
    'CUC': 'Leu', 'CCC': 'Pro', 'CAC': 'His', 'CGC': 'Arg',
    'CUA': 'Leu', 'CCA': 'Pro', 'CAA': 'Gln', 'CGA': 'Arg',
    'CUG': 'Leu', 'CCG': 'Pro', 'CAG': 'Gln', 'CGG': 'Arg',
    #A
    'AUU': 'Ile', 'ACU': 'Thr', 'AAU': 'Asn', 'AGU': 'Ser',
    'AUC': 'Ile', 'ACC': 'Thr', 'AAC': 'Asn', 'AGC': 'Ser',
    'AUA': 'Ile', 'ACA': 'Thr', 'AAA': 'Lys', 'AGA': 'Arg',
    'AUG': 'Met', 'ACG': 'Thr', 'AAG': 'Lys', 'AGG': 'Arg',
    #G
    'GUU': 'Val', 'GCU': 'Ala', 'GAU': 'Asp', 'GGU': 'Gly',
    'GUC': 'Val', 'GCC': 'Ala', 'GAC': 'Asp', 'GGC': 'Gly',
    'GUA': 'Val', 'GCA': 'Ala', 'GAA': 'Glu', 'GGA': 'Gly',
    'GUG': 'Val', 'GCG': 'Ala', 'GAG': 'Glu', 'GGG': 'Gly'
    }
    dnaCodonTable = {key.replace('U','T'):value for key, value in rnaCodonTable.items()}

def main():

    t = Converter()
        #take input and convert it to upper case
    code = (input("Please enter a codon or Amino Acid 1 or 3 letter code: ")).upper()

        #check if the input is in any of the dictionaries
    if code in t.short_AA:        #searches one letter dic
        print ('{} = {}'.format(code, t.short_AA[code].upper()))
    elif code in t.long_AA:         #searches three letter dic
        print ('{} = {}'.format(code, t.long_AA[code].upper()))
    elif code in t.rnaCodonTable:  #searches RNA codons)
        print ('{} = {}'.format(code, t.rnaCodonTable[code].upper()))
    elif code in t.dnaCodonTable:     #searches DNA codon table
        print ('{} = {}'.format(code, t.dnaCodonTable[code].upper()))
    else:
        print('unknown')        #if code is not in the dictionaries, print unknown
main()

