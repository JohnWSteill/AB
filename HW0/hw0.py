# John Steill 


import argparse, os, sys
import numpy as np


def read_file(input_filename):
    with open(input_filename, "r") as in_file:
        # Read the file contents
        # For normal NumPy data, you can load by np.load.
        # But for practice, you need to read them from plain text here.
        def is_header(line):
            return line[0] >= 'A' and line[0] <= 'Z'
        def is_end(line):
            return not line
        def text_to_array(line):
            return np.asarray([float(el) for el in line.split()])
        def arrays_to_matrix(arrays):
            return np.vstack(arrays)

        line = in_file.readline()
        matrices = {}
        while line:
            assert is_header(line), "Expected a character to start new matrix"
            key = line.strip()
            line = in_file.readline().strip()
            arrays = []
            while line:
                arrays.append(text_to_array(line))
                line = in_file.readline().strip()
            matrices[key] = arrays_to_matrix(arrays)
            line = in_file.readline()

        return matrices

def hw0_test(input_file, output_file):
    # out is the output file you can write to with the print function.
    with open(output_file, "w") as out: 
        def printf(s):
            print(s,file=out)
        # Part 1: String 
        # Create variables first_name and last_name with your name 
        first_name, last_name = ("John", "Steill")
    
        # Create full_name by concatenating the variables with a tab separating
        # them and print it
        full_name = '\t'.join((first_name, last_name))
        printf(full_name)

        # Transform the full_name string to all upper case letters and print it    
        printf(full_name.upper())
        
        # Part 2: list
        k = 4
        # Initialize a list x with k zeroes.
        x = [0 for x in range(k)]
    
        # Add one 1 to the head and append one 1 to the tail and print.
        x = [1] + x + [1]
        printf(x)

        # Set y to be the first 4 elements of x.
        y = x[:4]
 
    
        # Change the second last element of y to be 2 and print y
        y[-2] = 2
        printf(y)
    
    
        # Write a function to calculate the product of the elements in y.
        # Pass (skip over) the element if it is 0. Print the result.
        def get_prod(z):
    	    prod = 1
    	    for val in z:
    		    if val:
    			    prod *= val
    	    return prod
        printf(get_prod(y))
    
    
        # Python strings can be indexed in the same manner as lists
        course_str = "Advanced Bioinformatics"
        
        # Find the index of the 'B' and print the substring containing 'B' and
        # the next two characters ('Bio')
        b_ind = course_str.find('B')
        printf(course_str[b_ind:b_ind+3])
    
        # Part 3: dictionary
        # Note: strings can use single or double quotes.  This is set syntax.
        # Create a dictionary called hash_map.
        # Map the char a-d to 1-4. Save the mapping in hash_map.
        off = ord('a') - 1
        hash_map = {el:ord(el)-off for el in "abcd"}
    
        # Check if "e" exists in the hash_map. If not, map it to 6.
        # Print all key-value pairs in format <key:value> like "a:1".
        if 'e' not in hash_map:
    	    hash_map['e'] = 6
        printf(hash_map)

        #for ky,v in hash_map.items():
    	#    printf("{}:{}".format(ky,v))
    
        # Change "e" to 5.
        hash_map['e'] = 5
        # Print all key-value pairs in format <key:value> like "a:1".
        printf(hash_map)
        #for ky,v in hash_map.items():
    	#    printf("{}:{}".format(ky,v))
    
        # Part 4: NumPy array
        # Create a k by k random matrix using NumPy and print it.
        r_mat = np.random.rand(k,k)
        printf('\n'.join([" ".join([str(e) for e in el]) for el in r_mat]))

    
        # Print the minimum value in the matrix.
        printf(r_mat.min())
    
    
        alpha = 3
        # Complete the read_file function and call it to
        # read matrices A, B, C, D from input_file into NumPy arrays.   
        mats = read_file(input_file)
        A = mats['A']
        B = mats['B']
        C = mats['C']
        D = mats['D']
        
    
        # Print the shape of A.
        printf(A.shape)
    
        # Calculate F = alpha*A*B and print it. Here, we suggest you
        # write a function print_matrix(matrix, output_file) to print.
        def print_matrix(matrix, output_file):
            print('\n'.join([" ".join([str(e) for e in el]) for el in matrix]),
                    file=out)
    
        # Calculate G = A - 2*alpha*C and print it.
        G = A - 2*alpha*C
        print_matrix(G,output_file)
    
        # Part 5: NumPy binning
        # Print the min, the max and the range of the first column of D.
        D_1 = D[:,0]
        out_str = "{} {} {:.3f}"
        printf(out_str.format(D_1.min(), D_1.max(), D_1.max()-D_1.min()))
    
        # Assign the elements in the first column to 5 bins with equal width
        # in increasing order. 
        # Calculate bin width by dividing the range by the number of bins.
        hist = np.histogram(D_1, bins=5)
    
        # Print the number of elements that fall into each bin.
        printf(hist[0])
    
        # Print the min, the max and the range of the last row of D.
        D_m1 = D[:,-1]
        printf(out_str.format(D_m1.min(), D_m1.max(), D_m1.max()-D_m1.min()))
    
        # Assign the values in the last row to 5 bins in increasing order 
        # such that each bin contains the same number of elements.
        D_m1.sort()
        n = int(len(D_m1)/5)
        bins = [D_m1[x:x+n] for x in range(0,len(D_m1),n)]
    
        # Print the smallest and the largest values in each bin.
        [printf("{} {}".format(el.min(), el.max())) for el in bins]

def main(args):

    input_file = args.inputfile
    output_file = args.outputfile

    hw0_test(input_file, output_file)


# Note: this syntax checks if the Python file is being run as the main program
# and will not execute if the module is imported into a different module
if __name__ == "__main__":
    # Note: this example shows named command line arguments.  See the argparse
    # documentation for positional arguments and other examples.
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    # Note: you can use ' or " for strings
    parser.add_argument('--inputfile',
                        help='Enter a input data file path.',
                        type=str,
                        default='')
    parser.add_argument('--outputfile',
                        help='Enter a output file path.',
                        type=str,
                        default='')

    args = parser.parse_args()
    # Note: this simply calls the main function above, which we could have
    # given any name
    main(args)
