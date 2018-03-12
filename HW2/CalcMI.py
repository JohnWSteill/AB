import argparse, os, sys
import numpy as np
from scipy import stats
import itertools as it
# If you write your classes in other Python files
# you can import them here.

#def get_discrete_expr
PSEUDO_COUNT = .1

def get_dens_bins(x,n):
    '''For the equal density binning strategy, you can use np.percentile to 
    determine the limits of each bin. np.digitize will then assign each value 
    to a bin. Set right=True so that a bin includes its right limit. If you 
    do exactly this then your output will match ours. If you try something 
    different that also attempts to balance the number of values across the 
    bins, your result may be slightly off, but I think that is acceptable.'''

    return ( [-9999999] + [np.percentile(x,i*100/n) for i in range(1,n+1)])

def main(args):
    # Parse input arguments 
    data_inp_file = args.dataset
    bin_num = args.bin_num 
    output_file_path = args.out 
    bin_str = args.str.lower() 

    # ignore headers (genes are unnamed) skiprows=1 
    # not using time, so first column ommitted. 
    expr = np.loadtxt(args.dataset, skiprows=1)[:,1:]
    n_genes = np.shape(expr)[1]

    MI = np.zeros((n_genes, n_genes))
    outstrings = []
    if bin_str == 'kernel':
        for gene_a, gene_b in it.combinations(range(n_genes),2):
            a =  expr[:,gene_a]
            b = expr[:,gene_b]
            # make a mesh grid
            A,B = np.mgrid[-.1:1.1:100j, -.1:1.1:100j]
            # flatten to 2 x 10k
            positions = np.vstack([A.ravel(), B.ravel()])
            values = np.vstack([a, b])
            kernel = stats.gaussian_kde(values)
            K = np.reshape(kernel(positions).T, A.shape)
            K += .001;
            # Adapted from https://gist.github.com/GaelVaroquaux/ead9898bd3c973c40429
            K *= .012 * .012
            s1 = np.sum(K, axis=0).reshape((-1, K.shape[0]))
            s2 = np.sum(K, axis=1).reshape((K.shape[1], -1))
            MI = (np.sum(K * np.log2(K)) 
                    - np.sum(s1 * np.log2(s1))
                    - np.sum(s2 * np.log2(s2)))
            # Still not sure whats wrong with my two algorithms below
            '''
            g, p, dof, expected = stats.chi2_contingency(
                K, lambda_="log-likelihood")
            #print(0.5 * g / K.sum() / np.log(2))

            K /= K.sum()
            MI = 0
            for i in range(100):
                for j in range(100):
                    MI += K[i,j]*np.log2(K[i,j]/K[i,:].sum()/K[:,j].sum())
            '''
            outstring = '({},{})\t{:.3f}\n'.format(gene_a+1,gene_b+1, MI)
            outstrings.append(outstring)
        with open(output_file_path,'w') as f:
            f.writelines(sorted(outstrings,
                key = lambda s: -float(s.split('\t')[1])))
        return
    for gene_a, gene_b in it.combinations(range(n_genes),2):
        if bin_str=='uniform': 
            discrete_expr = np.histogram2d(
                    expr[:,gene_a],expr[:,gene_b],bins=bin_num)[0]
        elif bin_str=='density':
            a =  expr[:,gene_a]
            b = expr[:,gene_b]
            a_bin =  get_dens_bins(expr[:,gene_a],bin_num)
            b_bin =  get_dens_bins(expr[:,gene_b],bin_num)

            discrete_expr = np.histogram2d(a,b,bins=(a_bin,b_bin))[0]

        discrete_expr += PSEUDO_COUNT
        g, p, dof, expected = stats.chi2_contingency(
                discrete_expr, lambda_="log-likelihood")
        MI[gene_a,gene_b] = 0.5 * g / discrete_expr.sum() / np.log(2)
        outstring = '({},{})\t{:.3f}\n'.format(gene_a+1,gene_b+1, MI[gene_a,gene_b])
        outstrings.append(outstring)

    with open(output_file_path,'w') as f:
        f.writelines(sorted(outstrings,
            key = lambda s: -float(s.split('\t')[1])))








	# Where you run your code.

# Note: this syntax checks if the Python file is being run as the main program
# and will not execute if the module is imported into a different module
if __name__ == "__main__":
    # Note: this example shows named command line arguments.  See the argparse
    # documentation for positional arguments and other examples.
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--dataset',
                        help='input gene expression data file path',
                        type=str,
                        default='')
    parser.add_argument('--bin_num',
                        help='number of bins, not used with kernel binning',
                        type=int,
                        default=2)
    parser.add_argument('--out',
                        help='MI output file path',
                        type=str,
                        default='')
    parser.add_argument('--str',
                        help='binning strategy',
                        type=str,
                        choices={'uniform', 'density', 'kernel'},
                        default='uniform')

    args = parser.parse_args()

    main(args)
