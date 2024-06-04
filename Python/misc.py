from xml.dom.minicompat import NodeList
import pandas as pd
import numpy as np
import time
import os
from dataClass import inputParams, solveParams, hurricaneData, networkData
import pickle

#miscellaneous functions 

#sample a markovian state
def MC_sample(current_state,hurricaneDataSet):
    K = hurricaneDataSet.K;
    states = np.arange(K)
    weights = hurricaneDataSet.P_joint[current_state-1, :]
    kk = np.random.choice(states, p=weights)
    return kk

#training termination check
def termination_check(iter, relative_gap, LB, start, cutviol_iter, solveParams):
    flag = 0
    Elapsed = time.time() - start
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
            relative_gap = (LB[iter-1]-LB[iter-1-solveParams.stall])/max(1e-10,abs(LB[iter-1-solveParams.stall]))
            if relative_gap < solveParams.cutviol:
                flag = 4
                print("the LB is not making significant progress")
    return flag, Elapsed

# function to create an out-of-sample for given initial state k_init
def create_OSpaths(k_init,hurricaneDataSet):
    osfname = f"./data/case-study/SC-network/deterministic/OOS{k_init}-D.csv"
    OS_paths = pd.read_csv(osfname, header=0).to_numpy()
    rr = OS_paths.shape[0]
    cc = OS_paths.shape[1]
    print("rr = ", rr);
    print("cc = ", cc);
    OS_paths = np.full((rr, cc), k_init)
    for s in range(rr):
        for t in range(1, cc):
            OS_paths[s][t] = MC_sample(OS_paths[s][t-1],hurricaneDataSet)+1
    df = pd.DataFrame(OS_paths)
    df.to_csv(osfname, index=False)


# function to create a set of in-sample sample paths for each starting state k
def create_ISpaths(hurricaneDataSets,sample_size,filename):
    # open a file, where you stored the pickled data
    # file = open('data/synthetic/in_sample_100.dat', 'rb')
    # in_sample_dict = pickle.load(file)
    # file.close()
    ISpaths = {};
    K = hurricaneDataSets.K;
    for k in range(K):
        if k not in hurricaneDataSets.absorbing_states:
            in_sample = [];
            for s in range(sample_size):
                sample = [];
                sample.append(k);
                k_init = k
                while k_init not in hurricaneDataSets.absorbing_states:
                    k_next = MC_sample(k_init+1,hurricaneDataSets);
                    sample.append(k_next);
                    k_init = k_next;
                in_sample.append(sample);
            ISpaths[k] = in_sample;
    dbfile = open(filename, 'ab')
    # source, destination
    pickle.dump(ISpaths, dbfile)                    
    dbfile.close()

#function to save lp files
def save_lp(lp, name):
    with open(f"./output/lp/{name}.lp", "w") as f:
        f.write(lp)

#function save an n×m matrix as an csv file
def save_csv(m_data, m_colname, f_dir, m_fname):
    fname = os.path.join(f_dir, f"{m_fname}.csv")
    df = pd.DataFrame(m_data, columns=m_colname)
    df.to_csv(fname, index=False)

