import numpy as np
import itertools

ch_to_n = {el:i for i,el in enumerate('ACGT')} 

def get_pwm(seqs,z1,W): 
    Nc = {el:0 for el in 'ACGT'} 
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

def get_log_prob(seqs,pwm,z1,W):
    log_pwm = np.log10(pwm)
    prob = 0
    L = len(seqs[0])
    for i,seq in enumerate(seqs):
        motif_pos = range(W)
        non_motif_range = itertools.chain(range(z1[i]),range(z1[i]+W,L))
        for j in non_motif_range:
            prob += log_pwm[ch_to_n[seq[j]],0]
        for j in motif_pos:
            prob += log_pwm[ch_to_n[seq[j+z1[i]]],j+1]
    return prob



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
p4 =  get_log_prob(seqs, get_pwm(seqs,z1,4), z1, 4)
p5 =  get_log_prob(seqs, get_pwm(seqs,z1,5), z1, 5)
print("The log likelihood for W=4 solution is {:.2f}".format(p4))
print("The log likelihood for W=5 solution is {:.2f}".format(p5))
print("The difference is {:.3f}".format(p5-p4))




