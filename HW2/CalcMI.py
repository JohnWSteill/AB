import argparse, os, sys
import numpy as np
from scipy import stats
import itertools as it
# If you write your classes in other Python files
# you can import them here.

#def get_discrete_expr
PSEUDO_COUNT = .1

def get_dens_bins(x,n):
    '''For the equal density binning strategy, you can use np.percentile to determine the limits of each bin. np.digitize will then assign each value to a bin. Set right=True so that a bin includes its right limit. If you do exactly this then your output will match ours. If you try something different that also attempts to balance the number of values across the bins, your result may be slightly off, but I think that is acceptable.'''

    k = np.argsort(x)
    eps = 1e-5
    bins = [x[k[0]]-eps]
    n_in_bin = len(x)/n
    i = n_in_bin
    while i < len(x)-1-eps:
        bins.append(.5*(x[k[int(i)]] + x[k[int(i)+1]]))
        i += n_in_bin
    bins.append(x[k[-1]] + eps)
    return bins



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
        print('({},{})\t{:.3f}'.format(gene_a+1,gene_b+1, MI[gene_a,gene_b]))







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
