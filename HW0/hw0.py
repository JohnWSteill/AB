# Note: this is the Python syntax for importing modules
# We sometimes use package as a synonym for module, and formally
# a package is a module with a path on the file system.
# Modules may have submodules, which are delimited with .
# i.e. module.submodule
import argparse, os, sys

# Note: 'as' provides an alias for a module so that you can
# access its functions with syntax like np.mean() instead of numpy.mean()
import numpy as np


# Note: this is the syntax for Python functions.  You do not need to
# specify a return type in advance.  Multiple values can be returned
# using tuples such as (value1, value2, value3).
def read_file(input_filename):
    # Note: informally, 'with' statements are used for file reading and
    # writing, among other things.  They guarantee that the file is properly
    # closed regardless of whether the code block runs normally or throws
    # an exception.  The line below opens the file in read mode and creates
    # a file object.
    with open(input_filename, "r") as in_file:
        # Read the file contents
        # For normal NumPy data, you can load by np.load.
        # But for practice, you need to read them from plain text here.

    return

def hw0_test(input_file, output_file):
    # out is the output file you can write to with the print function.
        with open(output_file, "w") as out:

    
        # Part 1: String
        # Create variables first_name and last_name with your name
	first_name, last_name = ("John", "Steill")
    
        # Create full_name by concatenating the variables with a tab separating
        # them and print it
	full_name = '\t'.join(first_name, last_name)
	out.write(full_name + '\n')


    
        # Transform the full_name string to all upper case letters and print it    
	out.write(full_name.upper() + '\n')
        
        # Part 2: list
        k = 4
        # Initialize a list x with k zeroes.
	x = [0 for x in range(k)]
    
        # Add one 1 to the head and append one 1 to the tail and print.
        x = [1] + x + [1]
	out.write(x)

        # Set y to be the first 4 elements of x.
	y = x[:4]
 
    
        # Change the second last element of y to be 2 and print y
	y[-2] = 2
	out.write(y)
    
    
        # Write a function to calculate the product of the elements in y.
        # Pass (skip over) the element if it is 0. Print the result.
        def get_prod(z):
		prod = 1
		for val in z:
			if val:
				prod *= val
		return val
	out.write(get_prod(y))
    
    
        # Python strings can be indexed in the same manner as lists
        course_str = "Advanced Bioinformatics"
        
        # Find the index of the 'B' and print the substring containing 'B' and
        # the next two characters ('Bio')
	b_ind = course_str.find('B')
	out.write(course_str[b_ind:b_ind+3])
    
        # Part 3: dictionary
        # Note: strings can use single or double quotes.  This is set syntax.
        keys = {"a", "b", "c", "d"}
        # Create a dictionary called hash_map.
        # Map the char a-d to 1-4. Save the mapping in hash_map.
	hash_map = {el:ord(el) for el in "abcd"}
    
        # Check if "e" exists in the hash_map. If not, map it to 6.
        # Print all key-value pairs in format <key:value> like "a:1".
	if 'e' not in hash_map:
		hash_map['e'] = 6
	for k,v in hash_map.items():
		out.write("{}:{}\n".format(k,v))
    
        # Change "e" to 5.
	hash_map['e'] = 5
        # Print all key-value pairs in format <key:value> like "a:1".
	for k,v in hash_map.items():
		out.write("{}:{}\n".format(k,v))
    
        # Part 4: NumPy array
        # Create a k by k random matrix using NumPy and print it.
        r_mat = np.rand(k,k)
    
        # Print the minimum value in the matrix.
	
    
        alpha = 3
        # Complete the read_file function and call it to
        # read matrices A, B, C, D from input_file into NumPy arrays.   
    
        # Print the shape of A.
    
        # Calculate F = alpha*A*B and print it. Here, we suggest you
        # write a function print_matrix(matrix, output_file) to print.
    
        # Calculate G = A - 2*alpha*C and print it.
    
        # Part 5: NumPy binning
        # Print the min, the max and the range of the first column of D.
    
        # Assign the elements in the first column to 5 bins with equal width
        # in increasing order. 
        # Calculate bin width by dividing the range by the number of bins.
    
        # Print the number of elements that fall into each bin.
    
        # Print the min, the max and the range of the last row of D.
    
        # Assign the values in the last row to 5 bins in increasing order 
        # such that each bin contains the same number of elements.
    
        # Print the smallest and the largest values in each bin.

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
