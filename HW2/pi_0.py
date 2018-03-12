import numpy as np

A = np.array([
    [0.0, 20000, 0],
    [0.1, 15427, 0],
    [0.2, 12893, 0],
    [0.3, 11382, 0], 
    [0.4, 9834, 0],
    [0.5, 8466, 0],
    [0.6, 7030, 0],
    [0.7, 5259, 0],
    [0.8, 3484, 0],
    [0.9, 1714, 0]])

for row in A:
    row[2] = row[1] / (20000 * (1 - row[0]))
    print( "{:.1f} & {:d} & {:.2f} \\\\".format(row[0], int(row[1]), row[2]))


B = np.array([
    [0.000003, 0],
    [0.000007, 0],
    [0.000013, 0],
    [0.000024, 0],
    [0.000028, 0],
    [0.000033, 0],
    [0.000046, 0],
    [0.000055, 0],
    [0.000096, 0],
    [0.000099, 0]])

m = 20000
pi0 = .8

for i,row in enumerate(B):
    row[1] = pi0 * m * row[0] / (i+1)
print (B)

for i,row in enumerate(B):
    row[1] = min(B[i:,1])
    print( "{:d} & {:08.6f} & {:.3f} \\\\".format(int(i+1), row[0], row[1]))

print(B)
