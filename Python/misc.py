import pandas as pd
import numpy as np
import os
from dataClass import inputParams, solveParams, hurricaneData, networkData

#miscellaneous functions 

#sample a markovian state
def MC_sample(current_state,hurricaneDataSet):
    K = hurricaneDataSet.Na*hurricaneDataSet.Nb*hurricaneDataSet.T;
    states = np.arange(K)
    weights = hurricaneData.P_joint[current_state, :]
    kk = np.random.choice(states, p=weights)
    return kk

#training termination check
def termination_check(iter, relative_gap, LB, start, cutviol_iter, solveParams):
    flag = 0
    Elapsed = time() - start
    if iter > solveParams.max_iter:
        flag = 1
        print("max iteration is reached")
    elif Elapsed > solveParams.time_limit:
        flag = 2
        print("time limit is reached")
    elif cutviol_iter > solveParams.cutviol_maxiter:
        flag = 3
        print("cut violation is reached")
    else:
        if iter > solveParams.stall:
            relative_gap = (LB[iter]-LB[iter-solveParams.stall])/max(1e-10,abs(LB[iter-solveParams.stall]))
            if relative_gap < solveParams.cutviol:
                flag = 4
                print("the LB is not making significant progress")
return flag, Elapsed

# function to create an out-of-sample for given initial state k_init
def create_OSpaths(k_init):
    osfname = f"./data/OOS{k_init}.csv"
    OS_paths = pd.read_csv(osfname, header=0).to_numpy()
    rr = OS_paths.shape[0]
    cc = OS_paths.shape[1]
    OS_paths = np.full((rr, cc), k_init)
    for s in range(rr):
        for t in range(1, cc):
            OS_paths[s][t] = MC_sample(OS_paths[s][t-1])
    df = pd.DataFrame(OS_paths)
    df.to_csv(osfname, index=False)

#function to save lp files
def save_lp(lp, name):
    with open(f"./output/lp/{name}.lp", "w") as f:
        f.write(lp)

#function save an n×m matrix as an csv file
def save_csv(m_data, m_colname, f_dir, m_fname):
    fname = os.path.join(f_dir, f"{m_fname}.csv")
    df = pd.DataFrame(m_data, columns=m_colname)
    df.to_csv(fname, index=False)

