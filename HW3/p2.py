import numpy as np
cnt = np.asarray((620.,405.,5500.,40.,180.))
t_len = np.asarray((40.,270.,1000.,120.,100.))
amb = np.asarray((390., 4100.))

ec = (cnt/t_len)
f = ec/sum(ec)

rescue = np.asarray((
    cnt[0] + amb[0] * (f[0]/(f[0] + f[2])),
    cnt[1] + amb[1] * (f[1]/(f[1] + f[3])),
    cnt[2] + amb[0] * (f[2]/(f[0] + f[2])),
    cnt[3] + amb[1] * (f[3]/(f[1] + f[3])),
    cnt[4]))

ec_r = rescue/t_len
f_r = ec_r/sum(ec)
print rescue
f_r * 1000/f_r[-1]
