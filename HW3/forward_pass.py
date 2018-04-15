#!/usr/bin/env python
import os
import argparse
import numpy as np
from scipy import signal 
from scipy.stats import logistic as sigmoid

seq_len = 500
conv1_dim = (5,15)
conv2_dim = (5,15)
dense_dim = 0;
npools = 13;
max_pool_width = 35;

def one_hot(seq):
    T = {el:v for v,el in enumerate('ACGT')}
    seq_in_num = [T[el] for el in seq]
    return np.transpose(np.eye(4)[np.asarray(seq_in_num)])

def convolve2d(big,little):
    W = big.shape[1]
    w = little.shape[1]
    return np.asarray(
            [sum(sum(big[:,i:i+w]*little)) 
                for i in range(W-w+1)])

def do_fp(W,case,seq):
    seq_enc = one_hot(seq)
    nch = conv1_dim[0]
    width = conv1_dim[1]
    conv1_out = np.zeros((nch,seq_len - conv1_dim[1] + 1))
    for i in range(nch):
        ch = 'conv1_ch{}'.format(i+1)
        conv1_out[i,:] = convolve2d( seq_enc, W[ch].reshape(4,width))
        conv1_out[i,:] = np.maximum(conv1_out[i,:] + W['conv1_bias'][i],0)
    new_len = seq_len - width + 1;
    nch = conv2_dim[0]
    width = conv1_dim[1]
    conv2_out = np.zeros((nch,new_len - conv2_dim[1] + 1))
    for i in range(nch):
        ch = 'conv2_ch{}'.format(i+1)
        conv2_out[i,:] = convolve2d(conv1_out, W[ch].reshape(nch,width))
        conv2_out[i,:] = np.maximum(conv2_out[i,:] + W['conv2_bias'][i],0)
    conv2_out = conv2_out[:,:npools * max_pool_width]
    max_pool_out = np.hstack(
            np.apply_along_axis(max,1,el) for el in conv2_out.reshape(5,13,35))
    prob_in = np.dot(W['dense_wgts'], max_pool_out) + W['dense_bias']
    prob = 1./(1+np.exp(-prob_in));
    return prob


def read_input(weights):
    W = {}
    for w in os.listdir(weights):
        w_pth = os.path.join(weights,w)
        w_name = w.split('.')[0]
        W[w_name] = np.fromfile(w_pth, sep = ' ')
    return W

def main(args):
    W= read_input(args.weights)
    with open(args.out, 'w') as fout:
        for test_file in [args.positive_sequences, args.negative_sequences]:
            with open(test_file, 'r') as f:
                while True:
                    case = f.readline().strip()[1:]
                    if not case:
                        break
                    seq = f.readline().strip()
                    fout.write(str(do_fp(W,case,seq)[0]) + '\n')

if __name__ == "__main__":
    '''forward_pass.py --positive_sequences=<positive_sequences> 
    --negative_sequences=<negative_sequences> --weights=<weights> --out=<out>
    '''
    parser = argparse.ArgumentParser(
            description=__doc__, 
            formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--positive_sequences', type=str, default='')
    parser.add_argument('--negative_sequences', type=str, default='')
    parser.add_argument('--weights', type=str, default='')
    parser.add_argument('--out', type=str, default='') 
    args = parser.parse_args()
    main(args)


