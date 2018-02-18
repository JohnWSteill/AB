import argparse, os, sys
import numpy as np
import itertools

PI = .7

# You can choose to write classes in other python files
# and import them here.

def  E_step_get_Z(pwm, sequences):
    pass

def M_step_get_pwm(z,sequences):
    pass

def get_log_prob(pwm,z,sequences):
    pass

def load_sequences(f):
    pass


def get_init_motif(W,pi,load_sequences):
    init_log_prob = -999999999.
    background = .25 * np.ones((4,1));
    init_probs = { el: (1-pi)/3 * np.ones((4,1)) for el in 'ACGT'}
    for i,z in enumerate(init_probs): 
        init_probs[z][i] = pi
    def get_pwm_from_motif(motif):
        return np.hstack(((background,) + tuple(init_probs[z] for z in motif)))
    def get_k_mers_of_str(s,k):
        return list(s[i:i+k] for i in range(len(s)-k+1))
    def get_all_kmers(seqs,k):
        return (el for el in set(sum([get_k_mers_of_str(s,k) for s in seqs]))) 

    L = len(sequences[0])
    len_Z = L - W + 1
    Z = np.ones(n,len_Z) / len_Z
    for motif in get_all_kmers(sequences,W):
        pwm = get_pwm_from_motif(motif)
        z = E_step_get_Z(pwm, sequences)
        pwm = M_step_get_pwm(z,sequences)
        prob = get_log_prob(pwm,z,sequences)
        if prob > best_prob:
            best_prob = prob
            best_pwm = pwm
            best_z = z
        pwm = best_pwm
        z = best_z
     return (pwm, z)  

def write_output(pwm,z,sequences,args):
    pass


def main(args):
    model_file_path = args.M
    position_file_path = args.P
    subseq_file_path = args.SUBSEQ

    sequences = load_sequences(args.SEQ)
    (pwm, z) = get_init_motif(args.W, PI, sequences)
    
    pwm_old = pwm
    converged = False
    while not converged:
        pwm = get_pwm_from_motif(z, sequences)
        z = E_step_get_Z(pwm, sequences)
        converged = np.linalg.norm(pwm-pwm_old) < 1e-10
        pwm_old = pwm

    write_output(pwm, z, sequences, args)


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
