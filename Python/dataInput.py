import numpy as np;
import pandas as pd;
import sys;
import csv;
from dataClass import inputParams, solveParams, hurricaneData, networkData;

def hurricaneInput(intensityFile, locationFile, landfallFile, inputParams):
    dissipate_option = inputParams.dissipate_option;

    # probability distributions:
    P_intensity = pd.read_csv(intensityFile).values;  # intensity MC
    P_location = pd.read_csv(locationFile).values;  # location MC
    P_landfall = pd.read_csv(landfallFile).values;  # landfall MC

    Na = P_intensity.shape[0];  # intensity MC number of states
    Nb = P_location.shape[0];  # location MC number of states
    Nc = P_landfall.shape[0];  # landfall MC number of states
    T = Nc;  # define T_max

    Tmin = 2;  # for now we just hard code it

    K = Na * Nb * Nc;  # number of state in the joint MC
    P_joint = np.zeros((K, K));  # initialize the joint probability distribution MC
    S = [None] * K;  # list with elements [intensity, location]
    absorbing_states = [];  # list of absorbing states

    k1 = 0;  # counter for the number of states
    for k in range(1, Na+1):
        for l in range(1, Nb+1):
            for f in range(1, Nc+1):
                k1 += 1;
                k2 = 0;
                for n in range(1, Na+1):
                    for m in range(1, Nb+1):
                        for j in range(1, Nc+1):
                            k2 += 1;
                            P_joint[k1-1, k2-1] = P_intensity[k-1, n-1] * P_location[l-1, m-1] * P_landfall[f-1, j-1];
                S[k1-1] = [k, l, f];
                if dissipate_option:
                    if k == 1 or f == Nc:
                        absorbing_states.append(k1-1);
                else:
                    if f == Nc:
                        absorbing_states.append(k1-1);

                #if k == 2 and l == 2 and f == 1:
                #    print("k1 = ", k1);  57
    
                #if k == 4 and l == 2 and f == 1:
                #    print("k1 = ", k1); 155
    
                #if k == 6 and l == 2 and f == 1:
                #    print("k1 = ", k1); 253

    # normalize the probabilities
    P_temp = np.copy(P_joint);
    for k in range(K):
        for kk in range(K):
            P_joint[k, kk] = P_temp[k, kk] / np.sum(P_temp[k, :]);
    
    # Get the smallest transition probability that is nonzero to give us a correct threshold to filter out impossible transitions
    P_jointVec = np.copy(P_joint);
    nonzero_probs = P_jointVec[P_jointVec != 0];
    smallestTransProb = np.min(nonzero_probs) * 0.5;

    # Create a complete set of reachable nodes over time, starting from the initial state k_init
    k_init = inputParams.k_init;
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

    # Create a list of scenarios, along with the probability of occurrence, for each transient state node in the nodeList (set of reachable nodes from the initial state k_init)
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

    hurricaneDataSet = hurricaneData(P_intensity, P_location, P_landfall, Na, Nb, T, Tmin, P_joint, S, absorbing_states, smallestTransProb, nodeLists, nodeScenList, nodeScenWeights);
    return hurricaneDataSet;
    

def networkInput(Ni,Nj,costScalingFactor,netNodesFile,netParamsFile,hurricaneDataSet):
    nodes = pd.read_csv(netNodesFile);
    states = hurricaneDataSet.states;
    K = len(states);
    T = hurricaneDataSet.T;
    Na = hurricaneDataSet.Na;
    Nb = hurricaneDataSet.Nb;

    # List for the coordinates of the different supply points
    SP = [list(row) for row in nodes.iloc[:Ni, [0, 1]].values]

    # List for the coordinates of the different demand points
    DP = [list(row) for row in nodes.iloc[:Nj, [2, 3]].values]

    # Create an empty dictionary to store the data
    netParams = {}

    # Open the CSV file and read its contents
    with open(netParamsFile, mode='r') as file:
        csv_reader = csv.reader(file)
        
        # Iterate through each row in the CSV file
        for row in csv_reader:
            # The first element in each row is the key, and the rest are values
            key = row[0]
            values = row[1:]
            
            # Store the data in the dictionary
            netParams[key] = values

    # Now translate the csv data into parameters to use here:
    # 'MDC': x-coord, y-coord
    MDC = [int(netParams['MDC'][0]), int(netParams['MDC'][1])];
    # 'xbound": x_low, x_up
    xlow = int(netParams['x_bound'][0]);
    xup = int(netParams['x_bound'][1]);
    # 'other': [fuel, base, invCostRatio, penCostRatio, salvageCostRatio, dmax, cmax]
    fuel = float(netParams['other'][0]);
    base = float(netParams['other'][1]);
    invCostRatio = float(netParams['other'][2]);
    penCostRatio = float(netParams['other'][3]);
    salvageCostRatio = float(netParams['other'][4]);
    dMax = float(netParams['other'][5]);
    cMax = float(netParams['other'][6]);
    location_x = [None]*len(netParams['location_x']);
    for i in range(len(netParams['location_x'])):
        location_x[i] = int(netParams['location_x'][i]);
    location_y = [None]*len(netParams['location_y']);
    for i in range(len(netParams['location_y'])):
        location_y[i] = int(netParams['location_y'][i]);

    N0 = Ni + 1;

    # Unit cost of transporting/rerouting items from MDC/SP i to/between SP i'
    cb = np.empty((N0, Ni, T))
    for i in range(1, N0 + 1):
        for ii in range(1, Ni + 1):
            for t in range(1, T + 1):
                if i < N0:
                    cb[i - 1, ii - 1, t - 1] = (
                        fuel * np.linalg.norm(np.array(SP[i - 1]) - np.array(SP[ii - 1]), 2)
                        * (1 + costScalingFactor * (t - 1))
                    )
                else:
                    cb[i - 1, ii - 1, t - 1] = (
                        fuel
                        * np.linalg.norm(np.array(MDC) - np.array(SP[ii - 1]), 2)
                        * (1 + costScalingFactor * (t - 1))
                    )

    # Unit cost of transporting items from MDC/SP i to/between a demand point j
    ca = np.empty((N0, Nj, T))
    for i in range(1, N0 + 1):
        for j in range(1, Nj + 1):
            for t in range(1, T + 1):
                if i < N0:
                    ca[i - 1, j - 1, t - 1] = (
                        fuel
                        * np.linalg.norm(np.array(SP[i - 1]) - np.array(DP[j - 1]), 2)
                        * (1 + costScalingFactor * (t - 1))
                    )
                else:
                    ca[i - 1, j - 1, t - 1] = (
                        fuel
                        * np.linalg.norm(np.array(MDC) - np.array(DP[j - 1]), 2)
                        * (1 + costScalingFactor * (t - 1))
                    )

    cp = np.empty(T)
    ch = np.empty((Ni, T))
    for t in range(1, T + 1):
        cp[t - 1] = base * (1 + costScalingFactor * (t - 1))
        ch[:, t - 1] = np.full(Ni, 0.05 * base)

    p = penCostRatio * base
    q = salvageCostRatio * base
    x_cap = nodes.iloc[:Ni, 4].values * (Nj / Ni)
    x_0 = np.zeros(Ni)

    # Demand data
    SCEN = []

    for k in range(1, K + 1):
        scen = np.zeros(Nj)
        a = states[k - 1][0]
        l = states[k - 1][1]

        predicted = 0.5 * (location_x[l - 1] + location_y[l - 1])
        xx_coord = max(xlow, min(predicted, xup))
        landfall = [xx_coord, 0]

        for j in range(1, Nj + 1):
            c_j = np.linalg.norm(np.array(landfall) - np.array(DP[j - 1]), 2)
            if c_j <= cMax:
                scen[j - 1] = (
                    dMax * (1 - (c_j / cMax)) * (a - 1) ** 2 / ((Na - 1) ** 2)
                )
            else:
                scen[j - 1] = 0

        SCEN.append(scen)

    networkDataSet = networkData(Ni, Nj, SP, DP, fuel, cb, ca, ch, cp, p, q, dMax, x_cap, x_0, SCEN);
    return networkDataSet;