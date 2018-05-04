#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import argparse
import numpy as np

def get_pieces(s):
    for i in range(len(s)-1):
        yield s[:i+1]
        yield s[i+1:]
    yield s

def add_to_spectrum(sp,k,v):
    if k in sp:
        sp[k] = max(sp[k],v)
    else:
        sp[k] = v
    
def main(args):
    with open(args.mass_table, 'r') as f:
        mass = {el[0]: float(el[1]) for el in 
                [el_in.strip().split() for el_in in f.readlines()]}

    with open(args.peptides,'r') as pfile, open(args.out,'w') as ofile: 
        for line in pfile:
            peptide = line.strip()
            spectrum = {} 
            for piece in get_pieces(peptide): 
                m = int(np.round(sum([mass[el] for el in piece])))+1 
                add_to_spectrum(spectrum, m-1, 25) 
                add_to_spectrum(spectrum, m, 50) 
                add_to_spectrum(spectrum, m+1, 25) 
            ofile.write(peptide + '\n')
            for mz, h in sorted(spectrum.items()):
                ofile.write("{}\t{}\n".format(mz,h))
            ofile.write('\n')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
            description=__doc__, 
            formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--mass_table', type=str, default='mass_table.txt')
    parser.add_argument('--peptides', type=str, default='peptides.txt')
    parser.add_argument('--out', type=str, default='out_spectrum.txt') 
    args = parser.parse_args()
    main(args)


