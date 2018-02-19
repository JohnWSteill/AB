import argparse, os, sys
import numpy as np
import itertools
from scipy.special import logsumexp

PI = .7
CHAR_MP = {i:el for el,i in enumerate('ACGT')} 
L = 0
n = 0
W = 0


# You can choose to write classes in other python files
# and import them here.

def  E_step_get_Z(pwm, Z, sequences):
    global CHAR_MP, L, n, W
    pwm = np.log(pwm)
    for i,seq in enumerate(sequences):
        for j in range(L-W+1):
            motif_range = range(j,j+W)
            non_motif_range = itertools.chain(range(j),range(j+W,L))
            Z[i,j] = 0
            for k in non_motif_range:
                Z[i,j] += pwm[CHAR_MP[seq[k]],0]
            for k in motif_range:
                Z[i,j] += pwm[CHAR_MP[seq[k]],k-j+1]
    log_prob = sum(np.apply_along_axis(logsumexp,axis=1,arr=Z))
    z_nrm = np.exp(np.apply_along_axis(logsumexp,axis=1,arr=Z))
    Z = np.exp(Z) / z_nrm[:, np.newaxis]
    return (Z, log_prob)
                
def M_step_get_pwm(pwm, Z,sequences):
    global CHAR_MP, L, n, W
    pwm = np.ones((4,W+1))
    for i,seq in enumerate(sequences):
        for j in range(L-W+1):
            motif_range = range(j,j+W)
            non_motif_range = itertools.chain(range(j),range(j+W,L))
            for k in non_motif_range:
                pwm[CHAR_MP[seq[k]],0] += Z[i,j]
            for k in motif_range:
                pwm[CHAR_MP[seq[k]],k-j+1] += Z[i,j]
    pwm /= pwm.sum(axis=0)
    return pwm


def load_sequences(f):
    global CHAR_MP, L, n, W
    with open(f,'r') as seq_file:
        return [el.strip() for el in seq_file.readlines()]


def get_init_motif(sequences):
    global CHAR_MP, L, n, W
    best_prob = -999999999.
    background = .25 * np.ones((4,1));
    init_probs = { el: (1-PI)/3 * np.ones((4,1)) for el in 'ACGT'}
    for i,ch in enumerate(init_probs): 
        init_probs[ch][i] = PI
    def get_pwm_from_motif(motif):
        return np.hstack(((background,) 
            + tuple(init_probs[ch] for ch in motif)))
    def get_k_mers_of_str(s,k):
        return list(s[i:i+k] for i in range(len(s)-k+1))
    def get_all_kmers(seqs,k):
        return list(set(itertools.chain(
            *[get_k_mers_of_str(s,k) for s in seqs])))


    L = len(sequences[0])
    n = len(sequences)
    len_Z = L - W + 1
    Z = np.ones((n,len_Z)) / len_Z
    for motif in get_all_kmers(sequences,W):
        pwm = get_pwm_from_motif(motif)
        Z, prob = E_step_get_Z(pwm, Z, sequences)
        #pwm = M_step_get_pwm(pwm, Z,sequences)
        #Z, prob = E_step_get_Z(pwm, Z, sequences)
        if prob > best_prob:
            best_pwm = M_step_get_pwm(pwm, Z,sequences)
            best_prob = prob
            #print(best_prob, best_pwm)
            best_Z = Z
    Z = best_Z
    pwm = best_pwm
    return (pwm, Z)  

def write_output(pwm,Z,sequences,args):
    global CHAR_MP, L, n, W
    subseq_file_path = args.SUBSEQ

    with open(args.M,'w') as  model_file_path: 
        ''' A	0.227	0.049	0.831	0.844	0.231	0.028	0.209 '''
        model = []
        for  (ch, pwm_i) in zip('ACGT',pwm):
            model.append("{}\t".format(ch) +  
                    "\t".join(map("{:.3f}".format, pwm_i)) + '\t\n')
        model_file_path.writelines(model)

    pos_max = np.argmax(Z, axis=1)

    with open(args.P, 'w') as position_file_path:
        position_file_path.write('\n'.join([str(el) for el in pos_max]))

    with open(args.SUBSEQ,'w') as subseq_file_path:
        for p,seq in zip(pos_max, sequences):
            subseq_file_path.write(seq[p:p+W]+'\n')


def main(args):
    global CHAR_MP, L, n, W
    W = args.W

    sequences = load_sequences(args.SEQ)
    (pwm, Z) = get_init_motif(sequences)

    pwm_old = pwm
    converged = False
    while not converged:
        Z, prob = E_step_get_Z(pwm, Z, sequences)
        pwm = M_step_get_pwm(pwm, Z,sequences)
        converged = np.linalg.norm(pwm-pwm_old) < 1e-8
        #print ( np.linalg.norm(pwm-pwm_old))
        pwm_old = pwm

    write_output(pwm, Z, sequences, args)

# Note: this syntax checks if the Python file is being run as the main program
# and will not execute if the module is imported into a different module
if __name__ == "__main__":
    # Note: this example shows named command line arguments.  See the argparse
    # documentation for positional arguments and other examples.
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--SEQ',
                        help='sequences file path.',
                        type=str,
                        default='')
    parser.add_argument('--W',
                        help='length of the motif',
                        type=int,
                        default=6)
    parser.add_argument('--M',
                        help='model output file path.',
                        type=str,
                        default='')
    parser.add_argument('--P',
                        help='position output file path.',
                        type=str,
                        default='')
    parser.add_argument('--SUBSEQ',
                        help='subsequence output file path.',
                        type=str,
                        default='')

    args = parser.parse_args()
    # Note: this simply calls the main function above, which we could have
    # given any name
    # Parse input arguments
    # TODO It is generally good practice to validate the input arguments, e.g.,
    # verify that the input and output filenames are provided
    main(args)
