import numpy as np

def get_pwm(seqs,z1,W): 
    Nc = {el:0 for el in 'ACGT'} 
    ch_to_n = {el:i for i,el in enumerate('ACGT')} 
    for seq in seqs: 
        for ch in seq: 
            Nc[ch] += 1 
            
    pwm = np.zeros((4,W+1)) 
    for i, seq in enumerate(seqs): 
        for j in range(W): 
            pwm[ch_to_n[seq[j + z1[i]]],j+1] += 1 

    for ch in ch_to_n: 
        pwm[ch_to_n[ch],0] = Nc[ch] - sum(pwm[ch_to_n[ch],1:W+1])

    pwm = pwm/pwm.sum(axis=0)
    return (pwm)

def get_log_likelihood(seqs,z1,W,pwm):


seqs = [
        'CATGTGAA',
        'CAGCAGGG',
        'ACCTCTTC',
        'CAGACATG',
        'ACCTATCG',
        'GCGGCAGT',
        'GTGTAGTT',
        'CCAGGAAG',
        'ATGACCGG',
        'GGATAGTA']

z1 = [1,3,1,0,2,1,0,1,2,3]

print (get_pwm(seqs,z1,4))
pwm = get_pwm(seqs,z1,5)
for i, ch in enumerate('ACGT'):
    print (' & '.join(
            [ch] 
            + ["{:.3f}".format(el) for el in pwm[i,1:]] 
            + ["{:.3f}".format(pwm[i,0])]))



