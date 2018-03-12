import matplotlib
matplotlib.use('Agg')
import argparse, os, sys
import matplotlib.pyplot as plt
import numpy as np



def main(args):
    # Parse input arguments
    MI_file_path = args.MI
    gold_file_path = args.gold
    image_file_path = args.name + ".png"

    predicted_pairs = []        # predicted edges sorted by MI in decreasing order
    true_pairs = []             # true edges

    with open(MI_file_path, 'r') as f:
        lines = f.readlines()

    for line in lines:
        line = line.split()
        pair = line[0].split(',')
        predicted_pairs.append([pair[0][1:], pair[1][:-1]])
    
    with open(gold_file_path, 'r') as f:
        lines = f.readlines()
    
    for line in lines:
        line = line.split()
        true_pairs.append([line[0], line[1]])

    FPR = []
    TPR = []

    n = len(predicted_pairs)
    dy = 1/len(true_pairs)
    dx = 1/(n-len(true_pairs))
    (x,y) = (0,0)
    for pred in predicted_pairs:
        if pred in true_pairs:
            y += dy
        else:
            x += dx
        FPR.append(x)
        TPR.append(y)


    # TODO: Calculate the coordinates of points on ROC
    # store the x coordinates in FPR
    # store the y coordinates in TPR


    
    # Compute AUROC with the trapezoidal rule
    area = np.trapz(y=TPR, x=FPR)
    
    fig = plt.figure()
    plt.plot(FPR, TPR, '-')
    plt.xlabel('FPR')
    plt.ylabel('TPR')
    plt.title('ROC Curve (AUROC = {0:.3f})'.format(area))
    plt.axis([0, 1, 0, 1])
    
    # Increase the image resolution to 300 dots per inch
    fig.savefig(image_file_path, dpi=300)

# Note: this syntax checks if the Python file is being run as the main program
# and will not execute if the module is imported into a different module
if __name__ == "__main__":
    # Note: this example shows named command line arguments.  See the argparse
    # documentation for positional arguments and other examples.
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--MI',
                        help='MI file path.',
                        type=str,
                        default='')
    parser.add_argument('--gold',
                        help='gold standard network file path.',
                        type=str,
                        default='')
    parser.add_argument('--name',
                        help='output image file name without the file type.',
                        type=str,
                        default='')

    args = parser.parse_args()

    main(args)
