from generate import linear_gaussian
from idp import IBP
import numpy as np

if __name__ == '__main__':
    data = linear_gaussian()
    ibp = IBP(data["X"])
    print(ibp.compute_elbo())
    ibp.update_a()
    print("a updated")
    ibp.update_v()
    print("v_updated")
    ibp.update_z()
    print("z_updated")
    print(ibp.compute_elbo())



