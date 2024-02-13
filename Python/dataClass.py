#!/usr/bin/env python
# coding: utf-8

import numpy as np;
import pandas as pd;
import sys;
import csv;

class inputParams:
  def __init__(self,dissipate_option,absorbing_option,k_init,nbOS):
    self.dissipate_option = dissipate_option;
    self.absorbing_option = absorbing_option;
    self.k_init = k_init;
    self.nbOS = nbOS;

class solveParams:
  def __init__(self,max_iter,stall,cutviol_maxiter,time_limit,cutviol):
    self.max_iter = max_iter;
    self.stall = stall;
    self.cutviol_maxiter = cutviol_maxiter;
    self.time_limit = time_limit;
    self.cutviol = cutviol;

class hurricaneData:
  def createNodes(self, k_init):
    # create a list of MC states in each stage, starting from initial state k_init
    K = self.K;
    absorbing_states = self.absorbing_states;
    P_joint = self.P_joint;
    smallestTransProb = self.smallestTransProb;

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
    self.nodeLists = nodeLists;

    # Create a list of scenarios, along with the probability of occurrence, for each transient state node in the nodeList (set of reachable nodes from the initial state k_init)
    T = self.T;
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

    self.nodeScenList = nodeScenList; 
    self.nodeScenWeights = nodeScenWeights; 

  def input_from_Syn(self, intensityFile, locationFile, landfallFile, inputParams):
    # hurricane data class generator from synthetic instances
    dissipate_option = inputParams.dissipate_option;

    # probability distributions:
    P_intensity = pd.read_csv(intensityFile).values;  # intensity MC
    P_location = pd.read_csv(locationFile).values;  # location MC
    P_landfall = pd.read_csv(landfallFile).values;  # landfall MC

    Na = P_intensity.shape[0];  # intensity MC number of states
    Nb = P_location.shape[0];  # location MC number of states
    Nc = P_landfall.shape[0];  # landfall MC number of states
    T = Nc;  # define T_max

    K = Na * Nb * Nc;  # number of state in the joint MC
    P_joint = np.zeros((K, K));  # initialize the joint probability distribution MC
    S = [None] * K;  # list with elements [intensity, location]
    absorbing_states = [];  # list of absorbing states

    print("Total # of states K = ", K);
    print("Total # of intensity MC states = ", Na);
    print("Total # of location MC states = ", Nb);
    print("Total # of stages = ", Nc);

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

    # now store everything in the class
    self.Na = Na;
    self.K = K;
    self.T = T;
    self.P_joint = P_joint;
    self.states = S;
    self.absorbing_states = absorbing_states;
    self.smallestTransProb = smallestTransProb;

    # Create a complete set of reachable nodes over time, starting from the initial state k_init
    k_init = inputParams.k_init;
    self.createNodes(k_init); # List of MC states in each stage

  def input_from_Case(self, intensityFile, trackProbFile, trackErrorFile, landfallFile):
    # hurricane data class generator from synthetic instances
    # Default: trackProbFile = "data/case-study/mc_track_transition_prob_at_t"
    # Default: trackErrorFile = "data/case-study/mc_track_mean_error_at_t"

    # probability distributions:
    P_intensity = pd.read_csv(intensityFile).values;  # intensity MC
    P_landfall = pd.read_csv(landfallFile).values;  # landfall MC

    Na = P_intensity.shape[0];  # intensity MC number of states
    Nc = P_landfall.shape[0];  # landfall MC number of states
    T = Nc;  # define T_max

    # Now read in the track MC: note that each stage has their own track MC
    trackMatrices = [];
    trackStates = [];
    for t in range(T-1):
        temp_name = trackProbFile+str(t)+".csv";
        temp_MC = pd.read_csv(temp_name).values
        temp_name2 = trackErrorFile+str(t+1)+".csv";
        temp_states = pd.read_csv(temp_name2).values[:,1];

        trackMatrices.append(temp_MC);
        trackStates.append(temp_states);

        if temp_MC.shape[1] != len(temp_states):
            print("Error in reading track MCs!");
            exit(0);

    ########################################################################################
    # Data processing into a joint MC
    # First, count the total number of possible states
    K = 1;
    for t in range(T-1):
	    K += Na*len(trackStates[t]);

    print("Total # of states K = ", K);
    print("Total # of intensity MC states = ", Na);
    print("Total # of location MC states = [", end = "");
    for t in range(T-1):
        print(len(trackStates[t]), end = ""); 
        if t != T-2:
            print(",", end="");
        else:
            print("]");
    print("Total # of stages = ", Nc);
    
    P_joint = np.zeros((K, K));  # initialize the joint probability distribution MC
    S = [None] * K;  # list with elements [intensity, location] (actual state labels, starting from 1) 
    absorbing_states = [];  # list of absorbing states (their indices, starting from 0)

    k1 = 0;  # counter for the number of states
    S[0] = [1,1,1];
 
    for t in range(2,Nc+1):
        for l in range(1,len(trackStates[t-2])+1):
            for k in range(1,Na+1):
                k1 += 1;
                S[k1] = [k, l, t]
                if t == Nc:
                    absorbing_states.append(k1);
 
    for k in range(K):
        if S[k][2] == T:
            P_joint[k,k] = 1; # absorbing
        else:
            for kk in range(K):
                if S[kk][1] <= np.shape(trackMatrices[S[k][2]-1])[1]:
                    P_joint[k,kk] = P_intensity[S[k][0]-1,S[kk][0]-1]*trackMatrices[S[k][2]-1][S[k][1]-1,S[kk][1]-1]*P_landfall[S[k][2]-1,S[kk][2]-1]

    # normalize the probabilities
    P_temp = np.copy(P_joint);
    for k in range(K):
        for kk in range(K):
            P_joint[k, kk] = P_temp[k, kk] / np.sum(P_temp[k, :]);
    
    # Get the smallest transition probability that is nonzero to give us a correct threshold to filter out impossible transitions
    P_jointVec = np.copy(P_joint);
    nonzero_probs = P_jointVec[P_jointVec != 0];
    smallestTransProb = np.min(nonzero_probs) * 0.5;

    # now store everything in the class
    self.Na = Na;
    self.K = K;
    self.T = T;
    self.P_joint = P_joint;
    self.states = S;
    self.absorbing_states = absorbing_states;
    self.smallestTransProb = smallestTransProb;
    
    # Create a complete set of reachable nodes over time, starting from the initial state k_init = 1
    k_init = 1
    self.createNodes(k_init); # List of MC states in each stage

class networkData:
  def __init__(self, Ni, Nj):
    self.Ni = Ni;
    self.Nj = Nj;
    self.N0 = self.Ni+1;

  def input_from_Syn(self,costScalingFactor,netNodesFile,netParamsFile,hurricaneDataSet):
    # network data class generator from synthetic instances
    nodes = pd.read_csv(netNodesFile);
    states = hurricaneDataSet.states;
    K = len(states);
    T = hurricaneDataSet.T;
    Na = hurricaneDataSet.Na;

    # List for the coordinates of the different supply points
    SP = [list(row) for row in nodes.iloc[:self.Ni, [0, 1]].values]

    # List for the coordinates of the different demand points
    DP = [list(row) for row in nodes.iloc[:self.Nj, [2, 3]].values]

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

    # Unit cost of transporting/rerouting items from MDC/SP i to/between SP i'
    cb = np.empty((self.N0, self.Ni, T))
    for i in range(1, self.N0 + 1):
        for ii in range(1, self.Ni + 1):
            for t in range(1, T + 1):
                if i < self.N0:
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
    ca = np.empty((self.N0, self.Nj, T))
    for i in range(1, self.N0 + 1):
        for j in range(1, self.Nj + 1):
            for t in range(1, T + 1):
                if i < self.N0:
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
    ch = np.empty((self.Ni, T))
    for t in range(1, T + 1):
        cp[t - 1] = base * (1 + costScalingFactor * (t - 1))
        ch[:, t - 1] = np.full(self.Ni, invCostRatio * base)

    p = penCostRatio * base
    q = salvageCostRatio * base
    x_cap = nodes.iloc[:self.Ni, 4].values * (self.Nj / self.Ni)
    x_0 = np.zeros(self.Ni)

    # Demand data
    SCEN = []

    for k in range(1, K + 1):
        scen = np.zeros(self.Nj)
        a = states[k - 1][0]
        l = states[k - 1][1]

        predicted = 0.5 * (location_x[l - 1] + location_y[l - 1])
        xx_coord = max(xlow, min(predicted, xup))
        landfall = [xx_coord, 0]

        for j in range(1, self.Nj + 1):
            c_j = np.linalg.norm(np.array(landfall) - np.array(DP[j - 1]), 2)
            if c_j <= cMax:
                scen[j - 1] = (
                    dMax * (1 - (c_j / cMax)) * (a - 1) ** 2 / ((Na - 1) ** 2)
                )
            else:
                scen[j - 1] = 0

        SCEN.append(scen)

    # now store everything in the class
    self.fuel = fuel;
    self.cb = cb;
    self.ca = ca;
    self.ch = ch;
    self.cp = cp;
    self.p = p;
    self.q = q;
    self.x_cap = x_cap;
    self.x_0 = x_0;
    self.SCEN = SCEN;

  def input_from_Case(self,costScalingFactor,netFolderPath,netParamsFile,hurricaneDataSet):
    # network data class generator from synthetic instances
    d_JI = pd.read_csv(netFolderPath+"/d_JI.csv").values;
    d_II = pd.read_csv(netFolderPath+"/d_II.csv").values;
    d_KI = pd.read_csv(netFolderPath+"/d_KI.csv").values;
    d_KJ = pd.read_csv(netFolderPath+"/d_KJ.csv").values;
    d_SJ = pd.read_csv(netFolderPath+"/d_SJ.csv").values;
    x_cap = pd.read_csv(netFolderPath+"/x_cap.csv").values[0];
    max_D = pd.read_csv(netFolderPath+"/demand.csv").values;

    Nj = np.shape(d_JI)[0]; # number of DPs
    Ni = np.shape(d_JI)[1]; # number of SPs
    N0 = Ni + 1;

    states = hurricaneDataSet.states;
    K = hurricaneDataSet.K;
    T = hurricaneDataSet.T;
    Na = hurricaneDataSet.Na;

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
    # 'other': [fuel, base, invCostRatio, penCostRatio, salvageCostRatio, cmax]
    fuel = float(netParams['other'][0]);
    base = float(netParams['other'][1]);
    invCostRatio = float(netParams['other'][2]);
    penCostRatio = float(netParams['other'][3]);
    salvageCostRatio = float(netParams['other'][4]);
    cMax = float(netParams['other'][5]);
    

    # Unit cost of transporting/rerouting items from MDC/SP i to/between SP i'
    cb = np.empty((N0, Ni, T))
    for i in range(N0):
        for ii in range(Ni):
            for t in range(T):
                if i < N0-1:
                    cb[i, ii, t] = fuel * d_II[i,ii] * (1 + costScalingFactor * t)
                else:
                    cb[i, ii, t] = fuel * d_KI[0,ii] * (1 + costScalingFactor * t)

    # Unit cost of transporting items from MDC/SP i to/between a demand point j
    ca = np.empty((N0, Nj, T))
    for i in range(N0):
        for j in range(Nj):
            for t in range(T):
                if i < N0-1:
                    ca[i, j, t] = fuel * d_JI[j,i] * (1 + costScalingFactor * t)
                else:
                    ca[i, j, t] = fuel * d_KJ[0,j] * (1 + costScalingFactor * t)

    cp = np.empty(T)
    ch = np.empty((Ni, T))
    for t in range(1, T + 1):
        cp[t - 1] = base * (1 + costScalingFactor * (t - 1))
        ch[:, t - 1] = np.full(Ni, invCostRatio * base)

    p = penCostRatio * base
    q = salvageCostRatio * base
    x_0 = np.zeros(Ni)

    # Demand data
    SCEN = []

    for k in range(1, K + 1):
        scen = np.zeros(Nj)
        a = states[k - 1][0]
        l = states[k - 1][1]

        for j in range(1, Nj + 1):
            dLandfall = d_SJ[l-1,j-1]
            if dLandfall <= cMax:
                scen[j - 1] = max_D[0,j-1] * (1 - (dLandfall / cMax)) * (a - 1) ** 2 / ((Na - 1) ** 2)
            else:
                scen[j - 1] = 0

        SCEN.append(scen)

    # now store everything in the class
    self.fuel = fuel;
    self.cb = cb;
    self.ca = ca;
    self.ch = ch;
    self.cp = cp;
    self.p = p;
    self.q = q;
    self.x_cap = x_cap;
    self.x_0 = x_0;
    self.SCEN = SCEN;

    # Note: need to update Ni, Nj, N0 since this is from the datafile
    self.Ni = Ni;
    self.Nj = Nj;
    self.N0 = self.Ni+1;

#Ni: number of supply points (without MDC).
#Nj: number of demand points.
#N0: number of supply points including MDC.
#x_cap: the capacity of each supply point.
#cb: unit cost of transporting/rerouting relief items from the MDC/SP i to/between different SP i' 
#ch: unit cost for holding an item at SP i.
#cp: unit cost for relief item procurement.
#-----------------------------------------------------------------------------------------
#ca: unit cost of transporting relief items from the MDC/SP i to/between a demand point j.
#p: penalty cost for failing to serve the demand.
#q: salvage cost for each unit of overstock.
#-----------------------------------------------------------------------------------------
#Na: number of states in the intensity MC
#Nb: number of states in the locations MC
#K: number of states in the joint distrubutions between intesities and locations
#P_intesity: transition probability matrix of intensity MC
#P_location: transition probability matrix of location MC
#P_joint: transition probability matrix of the joint distrubution beween intensity and location MC
