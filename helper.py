import numpy as np

def arg_signchange(a):
    a_sign = np.sign(a)
    if_signflip = ((np.roll(a_sign, 1) - a_sign) != 0).astype(int)
    if_signflip[0] = 0
    arg_signflip = np.where(if_signflip == 1)
    return arg_signflip
