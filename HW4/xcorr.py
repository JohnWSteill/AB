#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import argparse
import numpy as np

def to_dense(s):
    s_d = np.zeros(max(s[:,0])+1)
    for (f,v) in s:
        s_d[int(f)] = v
    return s_d

def get_xcorr(exp, thr):
    n_exp, n_thr = (len(exp), len(thr))
    mx = max((n_exp, n_thr))
    exp_zp = np.zeros(mx+22)
    thr_zp = np.zeros(mx+22)
    exp_zp[10:10+n_exp] = exp
    thr_zp[10:10+n_thr] = thr

    R0 = np.dot(exp_zp,thr_zp)
    Rt = 0
    '''
    for j in range(-10,0):
        Rt += np.dot(thr_zp[-j:], exp_zp[:j])
    for j in range(1,11):
        Rt += np.dot(thr_zp[:-j], exp_zp[j:])
    return R0 - Rt/21.
    '''
    
    for i in range(10,mx-11):
        for j in range(-10,11):
            Rt += thr_zp[i]*exp_zp[i+j]
    return R0 - Rt/21.

def get_spectrum_from_library(fn):
    with open (fn,'r') as f:
        library = ''.join(f.readlines())
    libraries = library.split('\n\n')
    for lib in libraries:
        if not lib:
            continue
        lib_c = lib.split('\n')
        peptide = lib_c[0]
        spectrum = to_dense(np.fromstring('\t'.join(lib_c[1:]),sep='\t').reshape(-1,2))
        yield (peptide, spectrum)

def main(args):
    query_spectrum = to_dense(np.fromfile(args.query, sep='\t').reshape(-1,2))
    all_corrs = [(get_xcorr(query_spectrum, spectrum), peptide) 
            for (peptide, spectrum) in get_spectrum_from_library(args.library)]
    all_corrs.sort(key = lambda x: -x[0])
    nrm = all_corrs[0][0]
    print "nrm: {}".format(nrm)
    all_corrs = [(el[0]/nrm ,el[1]) for el in all_corrs]
    with open(args.out, 'w') as f:
        f.writelines(["{:.3f}\t{}\n".format(el[0],el[1]) for el in all_corrs])

if __name__ == "__main__":
    ''' xcorr.py --query=<spectrum> --library=<library> --out=<out> '''
    parser = argparse.ArgumentParser(
            description=__doc__, 
            formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--query', type=str, default='spectrum1.txt')
    parser.add_argument('--library', type=str, default='out_spectrum.txt')
    parser.add_argument('--out', type=str, default='out_xcorr1.txt') 
    args = parser.parse_args()
    main(args)


