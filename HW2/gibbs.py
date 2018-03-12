import numpy as np

adj_mat = np.zeros((15,15))
adj_mat[0,1] = adj_mat[1,0] = adj_mat[0,2] = adj_mat[2,0] = 1
adj_mat[1,3] = adj_mat[3,1] = adj_mat[1,4] = adj_mat[4,1] = 1
adj_mat[2,5] = adj_mat[5,2] = adj_mat[2,6] = adj_mat[6,2] = 1
adj_mat[3,7] = adj_mat[7,3] = adj_mat[3,8] = adj_mat[8,3] = 1
adj_mat[4,9] = adj_mat[9,4] = adj_mat[4,10] = adj_mat[10,4] = 1
adj_mat[5,11] = adj_mat[11,5] = adj_mat[5,12] = adj_mat[12,5] = 1
adj_mat[6,13] = adj_mat[13,6] = adj_mat[6,14] = adj_mat[14,6] = 1

State = -1 * np.ones(15)
State[14] = 1;

match_prob = np.exp(1)
nomatch_prob = np.exp(-1)
np.seterr(all='raise')
def gibbs_iter(state):
    node = np.random.random_integers(0,13)
    pairs = adj_mat[node,:]*state
    on = np.count_nonzero(pairs==1)
    off = np.count_nonzero(pairs==-1)
    p_on = match_prob * on + nomatch_prob * off
    p_off = match_prob *off + nomatch_prob * on
    p_tot = p_on + p_off
    gibbs_magic = np.random.uniform()
    if gibbs_magic < p_off/p_tot:
        state[node] = -1
    else:
        state[node] = 1
    return (state)

for i in range(10000):
    State = gibbs_iter(State) 
count = 0
for i in range(50000):
    State = gibbs_iter(State) 
    if State[0]==1 :
        count += 1
print (count/50000)



