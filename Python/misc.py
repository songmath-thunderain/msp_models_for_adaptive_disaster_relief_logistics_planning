from xml.dom.minicompat import NodeList
import pandas as pd
import numpy as np
import time
import os
from dataClass import inputParams, solveParams, hurricaneData, networkData

#miscellaneous functions 

#sample a markovian state
def MC_sample(current_state,hurricaneDataSet):
    K = hurricaneDataSet.K;
    states = np.arange(K)
    weights = hurricaneDataSet.P_joint[current_state, :]
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
    osfname = f"./data/case-study/OOS{k_init}.csv"
    OS_paths = pd.read_csv(osfname, header=0).to_numpy()
    rr = OS_paths.shape[0]
    cc = OS_paths.shape[1]
    OS_paths = np.full((rr, cc), k_init)
    for s in range(rr):
        for t in range(1, cc):
            OS_paths[s][t] = MC_sample(OS_paths[s][t-1],hurricaneDataSet)
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

# create a list of MC states in each stage, starting from initial state k_init
def createNodes(k_init, K, absorbing_states, P_joint, smallestTransProb):
    nodeLists = []
    nodeLists.append([k_init-1])
    stopFlag = False

    while not stopFlag:
        tempList = []
        stopFlag = True

        for k in range(K):
            for kk in nodeLists[-1]:
                if (kk not in absorbing_states) and (P_joint[kk][k] > smallestTransProb):
                    tempList.append(k)

                    if k not in absorbing_states:
                        stopFlag = False
                    break

        nodeLists.append(tempList)
    return nodeLists

# Create a list of scenarios, along with the probability of occurrence, for each transient state node in the nodeList (set of reachable nodes from the initial state k_init)
def createNodeScens(nodeLists, T, absorbing_states, P_joint, smallestTransProb):
    nodeScenList = {}
    nodeScenWeights = {}

    for t in range(T - 2, -1, -1):
        # Starting from T-2 since at T-1, all states should be absorbing
        for k in range(len(nodeLists[t])):
            if nodeLists[t][k] not in absorbing_states:
                nodeScenList[(t, nodeLists[t][k])] = []
                nodeScenWeights[(t, nodeLists[t][k])] = []

                for kk in range(len(nodeLists[t + 1])):
                    if P_joint[nodeLists[t][k]][nodeLists[t + 1][kk]] > smallestTransProb:
                        if nodeLists[t + 1][kk] in absorbing_states:
                            # absorbing states, directly append
                            nodeScenList[(t, nodeLists[t][k])].append((t + 1, nodeLists[t + 1][kk]))
                            nodeScenWeights[(t, nodeLists[t][k])].append(P_joint[nodeLists[t][k]][nodeLists[t + 1][kk]])
                        else:
                            # transient states, append the corresponding scenlist and weights
                            for j in range(len(nodeScenList[(t + 1, nodeLists[t + 1][kk])])):
                                if (nodeScenList[(t + 1, nodeLists[t + 1][kk])][j] not in nodeScenList[(t, nodeLists[t][k])]):
                                    # Not in the scenario list, so go ahead and add it
                                    nodeScenList[(t, nodeLists[t][k])].append(nodeScenList[(t + 1, nodeLists[t + 1][kk])][j])
                                    nodeScenWeights[(t, nodeLists[t][k])].append(P_joint[nodeLists[t][k]][nodeLists[t + 1][kk]] * nodeScenWeights[(t + 1, nodeLists[t + 1][kk])][j])
                                else:
                                    # in the scenario list, increment the probability
                                    ind = list(nodeScenList[(t, nodeLists[t][k])]).index(next((x for x in nodeScenList[(t, nodeLists[t][k])] if x == nodeScenList[(t + 1, nodeLists[t + 1][kk])][j]), None))
                                    nodeScenWeights[(t, nodeLists[t][k])][ind] += P_joint[nodeLists[t][k]][nodeLists[t + 1][kk]]*nodeScenWeights[(t + 1, nodeLists[t + 1][kk])][j]

                if abs(sum(nodeScenWeights[(t, nodeLists[t][k])]) - 1) > 1e-6:
                    print("Wrong!")
                    sys.exit(0)

    return nodeScenList, nodeScenWeights