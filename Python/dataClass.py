#!/usr/bin/env python
# coding: utf-8

import numpy as np;
import pandas as pd;
import sys;
import csv;
import json;
import ast;
from shapely.geometry import Point, LineString

class inputParams:
  def __init__(self,dissipate_option,absorbing_option,cost_structure,safe_time,tau,k_init,nbOS):
    self.dissipate_option = dissipate_option;
    self.absorbing_option = absorbing_option;
    self.cost_structure = cost_structure;
    self.safe_time = safe_time;
    self.k_init = k_init;
    self.tau = tau;
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
    # create a complete list of reachable node over time, starting from initial state k_init
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
        # if there is at least one state that is not absorbing, turn the stopFlag back to false
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
    # Note: there might be repetitions in the nodeLists, e.g., the same Hurricane state can be reached at different times
    T = self.T;
    nodeScenList = {}
    nodeScenWeights = {}
    nodeTime2Go = {}

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

                time2Go = 0;
                for l in range(len(nodeScenList[(t, nodeLists[t][k])])):
                    time2Go += (nodeScenList[(t, nodeLists[t][k])][l][0]-t)*nodeScenWeights[(t, nodeLists[t][k])][l];
                nodeTime2Go[(t, nodeLists[t][k])] = time2Go;
    self.nodeScenList = nodeScenList; 
    self.nodeScenWeights = nodeScenWeights; 
    self.nodeTime2Go = nodeTime2Go;

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
    print("Total # of absorbing stages = ", len(absorbing_states));

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
    # hurricane data class generator from case study instances
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
  
  def input_from_Case_new(self, absorbingFile, MCFile):
    # hurricane data class generator from case study instances
    with open(MCFile, 'r') as file:
        MC = json.load(file)

    if absorbingFile != None:
        with open(absorbingFile, 'r') as file:
            absorbing = json.load(file)
    T = len(list(MC.keys()));  # define T_max

    ########################################################################################
    # Data processing into a joint MC
    # First, count the total number of possible states
    K = 1;
    for t in range(T-1):
	    K += len(list(MC[str(t+1)].keys()));

    print("Total # of states K = ", K);
    print("Total # of MC states = [", end = "");
    for t in range(T-1):
        print(len(list(MC[str(t+1)].keys())), end = ""); 
        if t != T-2:
            print(",", end="");
        else:
            print("]");
    print("Total # of stages = ", T);
    
    P_joint = np.zeros((K, K));  # initialize the joint probability distribution MC
    S = [None] * K;  # list with elements [intensity, location] (actual state labels) 
    absorbing_states = [];  # list of absorbing states (their indices, starting from 0)

    k1 = 0;  # counter for the number of states
    initS = ast.literal_eval(list(MC['0'].keys())[0])
    S[0] = [initS[1]+1,initS[0],1]; # note that the intensity state starts from 1
 
    for t in range(2,T+1):
        for k in list(MC[str(t-1)].keys()):
            k1 += 1;
            tempS = ast.literal_eval(k);
            S[k1] = [tempS[1]+1, tempS[0], t]
            if absorbingFile == None:
                # absorbingFile is not provided, deterministic case where all states in t = T are absorbing
                if t == T:
                    absorbing_states.append(k1);
            else:
                # absorbingFile is provided, get the absorbing states accordingly
                if absorbing[str(t-1)][k]:
                    absorbing_states.append(k1);
 
    for k in range(K):
        if k in absorbing_states:
            P_joint[k,k] = 1; # absorbing
        else:
            for kk in range(K):
                if S[kk][2] == S[k][2] + 1:
                    P_joint[k,kk] = MC[str(S[k][2]-1)]['('+str(S[k][1])+', '+str(S[k][0]-1)+')']['('+str(S[kk][1])+', '+str(S[kk][0]-1)+')']
                else:
                    P_joint[k,kk] = 0;    

    '''
    # normalize the probabilities
    P_temp = np.copy(P_joint);
    for k in range(K):
        for kk in range(K):
            P_joint[k, kk] = P_temp[k, kk] / np.sum(P_temp[k, :]);
    '''
    # Get the smallest transition probability that is nonzero to give us a correct threshold to filter out impossible transitions
    P_jointVec = np.copy(P_joint);
    nonzero_probs = P_jointVec[P_jointVec != 0];
    smallestTransProb = np.min(nonzero_probs) * 0.5;

    # now store everything in the class
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

  def input_from_Syn(self,cost_structure,safe_time,costScalingFactor,netNodesFile,netParamsFile,hurricaneDataSet):
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
    cb = np.empty((self.N0, self.Ni, T, K))
    for i in range(1, self.N0 + 1):
        for ii in range(1, self.Ni + 1):
            for t in range(1, T + 1):
                if i < self.N0:
                    if cost_structure == 0:
                        # cost_structure is only time dependent
                        for k in range(K):
                            cb[i - 1, ii - 1, t - 1, k] = (
                            fuel * np.linalg.norm(np.array(SP[i - 1]) - np.array(SP[ii - 1]), 2)
                            * (1 + costScalingFactor * (t - 1))
                            )
                    if cost_structure == -1 or cost_structure == 1:
                        for k in range(K):
                            surgeFlag = False;
                            if (t-1,k) in hurricaneDataSet.nodeTime2Go:
                                if hurricaneDataSet.nodeTime2Go[(t-1,k)] <= safe_time + 1e-5:
                                    surgeFlag = True;
                            else:
                                surgeFlag = True;
                            if not surgeFlag: 
                                cb[i - 1, ii - 1, t - 1, k] = fuel * np.linalg.norm(np.array(SP[i - 1]) - np.array(SP[ii - 1]), 2)
                            else:
                                cb[i - 1, ii - 1, t - 1, k] = fuel * np.linalg.norm(np.array(SP[i - 1]) - np.array(SP[ii - 1]), 2) * costScalingFactor
                else:
                    if cost_structure == 0:
                        # cost_structure is only time dependent
                        for k in range(K):
                            cb[i - 1, ii - 1, t - 1, k] = (
                            fuel
                            * np.linalg.norm(np.array(MDC) - np.array(SP[ii - 1]), 2)
                            * (1 + costScalingFactor * (t - 1))
                            )
                    if cost_structure == -1 or cost_structure == 1:
                        for k in range(K):
                            surgeFlag = False;
                            if (t-1,k) in hurricaneDataSet.nodeTime2Go:
                                if hurricaneDataSet.nodeTime2Go[(t-1,k)] <= safe_time + 1e-5:
                                    surgeFlag = True;
                            else:
                                surgeFlag = True;
                            if not surgeFlag: 
                                cb[i - 1, ii - 1, t - 1, k] = (
                                fuel
                                * np.linalg.norm(np.array(MDC) - np.array(SP[ii - 1]), 2))
                            else:
                                cb[i - 1, ii - 1, t - 1, k] = (
                                fuel
                                * np.linalg.norm(np.array(MDC) - np.array(SP[ii - 1]), 2)
                                * costScalingFactor)
    # Unit cost of transporting items from MDC/SP i to/between a demand point j
    ca = np.empty((self.N0, self.Nj, T, K))
    for i in range(1, self.N0 + 1):
        for j in range(1, self.Nj + 1):
            for t in range(1, T + 1):
                if i < self.N0:
                    if cost_structure == 0:
                        # cost_structure is only time dependent
                        for k in range(K):
                            ca[i - 1, j - 1, t - 1, k] = (
                            fuel
                            * np.linalg.norm(np.array(SP[i - 1]) - np.array(DP[j - 1]), 2)
                            * (1 + costScalingFactor * (t - 1))
                            )
                    if cost_structure == -1 or cost_structure == 1:
                        for k in range(K):
                            surgeFlag = False;
                            if (t-1,k) in hurricaneDataSet.nodeTime2Go:
                                if hurricaneDataSet.nodeTime2Go[(t-1,k)] <= safe_time + 1e-5:
                                    surgeFlag = True;
                            else:
                                surgeFlag = True;
                            if not surgeFlag:
                                ca[i - 1, j - 1, t - 1, k] = fuel * np.linalg.norm(np.array(SP[i - 1]) - np.array(DP[j - 1]), 2)
                            else:
                                ca[i - 1, j - 1, t - 1, k] = fuel * np.linalg.norm(np.array(SP[i - 1]) - np.array(DP[j - 1]), 2)*costScalingFactor    
                else:
                    if cost_structure == 0:
                        # cost_structure is only time dependent
                        for k in range(K):
                            ca[i - 1, j - 1, t - 1, k] = (
                            fuel
                            * np.linalg.norm(np.array(MDC) - np.array(DP[j - 1]), 2)
                            * (1 + costScalingFactor * (t - 1))
                            )
                    if cost_structure == -1 or cost_structure == 1:
                        for k in range(K):
                            surgeFlag = False;
                            if (t-1,k) in hurricaneDataSet.nodeTime2Go:
                                if hurricaneDataSet.nodeTime2Go[(t-1,k)] <= safe_time + 1e-5:
                                    surgeFlag = True;
                            else:
                                surgeFlag = True;
                            if not surgeFlag:
                                ca[i - 1, j - 1, t - 1, k] = fuel * np.linalg.norm(np.array(MDC) - np.array(DP[j - 1]), 2)
                            else:
                                ca[i - 1, j - 1, t - 1, k] = fuel * np.linalg.norm(np.array(MDC) - np.array(DP[j - 1]), 2)*costScalingFactor
    cp = np.empty((T, K))
    ch = np.empty((self.Ni, T))
    for t in range(1, T + 1):
        if cost_structure == 0:
            # cost_structure is only time dependent
            for k in range(K):
                cp[t - 1, k] = base * (1 + costScalingFactor * (t - 1))
        if cost_structure == -1 or cost_structure == 1:
            for k in range(K):
                surgeFlag = False;
                if (t-1,k) in hurricaneDataSet.nodeTime2Go:
                    if hurricaneDataSet.nodeTime2Go[(t-1,k)] <= safe_time + 1e-5:
                        surgeFlag = True;
                else:
                    surgeFlag = True;
                if not surgeFlag:
                    cp[t-1, k] = base
                else:
                    cp[t-1, k] = base*costScalingFactor;
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

  def input_from_Case(self,cost_structure,safe_time,costScalingFactor,netFolderPath,netParamsFile,hurricaneDataSet):
    # data input interface for case study (old format)
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
    cb = np.empty((N0, Ni, T, K))
    for i in range(N0):
        for ii in range(Ni):
            for t in range(T):
                if i < N0-1:
                    if cost_structure == 0:
                        # cost_structure is only time dependent
                        for k in range(K):
                            cb[i, ii, t, k] = fuel * d_II[i,ii] * (1 + costScalingFactor * t)
                    if cost_structure == -1 or cost_structure == 1:
                        for k in range(K):
                            surgeFlag = False;
                            if (t,k) in hurricaneDataSet.nodeTime2Go:
                                if hurricaneDataSet.nodeTime2Go[(t,k)] <= safe_time + 1e-5:
                                    surgeFlag = True;
                            else:
                                surgeFlag = True;
                            if not surgeFlag:
                                cb[i, ii, t, k] = fuel * d_II[i,ii]
                            else:
                                cb[i, ii, t, k] = fuel * d_II[i,ii]*costScalingFactor
                else:
                    if cost_structure == 0:
                        # cost_structure is only time dependent
                        for k in range(K):
                            cb[i, ii, t, k] = fuel * d_KI[0,ii] * (1 + costScalingFactor * t)
                    if cost_structure == -1 or cost_structure == 1:
                        for k in range(K):
                            surgeFlag = False;
                            if (t,k) in hurricaneDataSet.nodeTime2Go:
                                if hurricaneDataSet.nodeTime2Go[(t,k)] <= safe_time + 1e-5:
                                    surgeFlag = True;
                            else:
                                surgeFlag = True;
                            if not surgeFlag:
                                cb[i, ii, t, k] = fuel * d_KI[0,ii]
                            else:
                                cb[i, ii, t, k] = fuel * d_KI[0,ii]*costScalingFactor
    # Unit cost of transporting items from MDC/SP i to/between a demand point j
    ca = np.empty((N0, Nj, T, K))
    for i in range(N0):
        for j in range(Nj):
            for t in range(T):
                if i < N0-1:
                    if cost_structure == 0:
                        # cost_structure is only time dependent
                        for k in range(K):
                            ca[i, j, t, k] = fuel * d_JI[j,i] * (1 + costScalingFactor * t)
                    if cost_structure == -1 or cost_structure == 1:
                        for k in range(K):
                            surgeFlag = False;
                            if (t,k) in hurricaneDataSet.nodeTime2Go:
                                if hurricaneDataSet.nodeTime2Go[(t,k)] <= safe_time + 1e-5:
                                    surgeFlag = True;
                            else:
                                surgeFlag = True;
                            if not surgeFlag:
                                ca[i, j, t, k] = fuel * d_JI[j,i]
                            else:
                                ca[i, j, t, k] = fuel * d_JI[j,i]*costScalingFactor
                else:
                    if cost_structure == 0:
                        # cost_structure is only time dependent
                        for k in range(K):
                            ca[i, j, t, k] = fuel * d_KJ[0,j] * (1 + costScalingFactor * t)
                    if cost_structure == -1 or cost_structure == 1:
                        for k in range(K):
                            surgeFlag = False;
                            if (t,k) in hurricaneDataSet.nodeTime2Go:
                                if hurricaneDataSet.nodeTime2Go[(t,k)] <= safe_time + 1e-5:
                                    surgeFlag = True;
                            else:
                                surgeFlag = True;
                            if not surgeFlag:
                                ca[i, j, t, k] = fuel * d_KJ[0,j]
                            else:
                                ca[i, j, t, k] = fuel * d_KJ[0,j]*costScalingFactor

    cp = np.empty((T, K))
    ch = np.empty((Ni, T))
    for t in range(1, T + 1):
        if cost_structure == 0:
            # cost_structure is only time dependent
            for k in range(K):
                cp[t - 1, k] = base * (1 + costScalingFactor * (t - 1))
        if cost_structure == -1 or cost_structure == 1:
            for k in range(K):
                surgeFlag = False;
                if (t-1,k) in hurricaneDataSet.nodeTime2Go:
                    if hurricaneDataSet.nodeTime2Go[(t-1,k)] <= safe_time + 1e-5:
                        surgeFlag = True;
                else:
                    surgeFlag = True;
                if not surgeFlag:
                    cp[t - 1, k] = base
                else:
                    cp[t - 1, k] = base*costScalingFactor
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

    print("# of SPs = ", self.Ni);
    print("# of DPs = ", self.Nj);

  def input_from_Case_new(self,cost_structure,safe_time,costScalingFactor,netFolderPath,netParamsFile,hurricaneDataSet):
    # data input interface for case study (new format)
    df = pd.read_excel(netFolderPath+'locations.xlsx');
    Nj = 0; # of DPs
    Ni = 0; # of SPs
    for i in range(df.shape[0]):
        if df['Type'][i] == 'RSA':
            Ni += 1;
        if df['Type'][i] == 'PoD':
            Nj += 1;
    N0 = Ni + 1;

    states = hurricaneDataSet.states;
    K = hurricaneDataSet.K;
    T = hurricaneDataSet.T;

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
    f_cap = float(netParams['other'][6]);

    d_II = {};
    d_JI = {};
    d_KI = {};
    d_KJ = {};

    # now propagate these dictionaries from the dataframe data and distance matrix
    # Read in the distance matrix
    distanceMatrix = np.loadtxt(netFolderPath+'distanceMatrix.txt');

    # WARNING: hardcode on the indices of the distance matrix
    for i in range(Ni):
        for ii in range(Ni):
            if i == ii:
                d_II[i,ii] = 0;
            else:
                d_II[i,ii] = distanceMatrix[i+1][ii+1];
    
    for j in range(Nj):
        for i in range(Ni):
            d_JI[j,i] = distanceMatrix[j+1+Ni][i+1];
    
    for i in range(Ni):
        d_KI[i] = distanceMatrix[0][i+1];
    
    for j in range(Nj):
        d_KJ[j] = distanceMatrix[0][j+1+Ni];
        
    # Unit cost of transporting/rerouting items from MDC/SP i to/between SP i'
    cb = np.empty((N0, Ni, T, K))
    for i in range(N0):
        for ii in range(Ni):
            for t in range(T):
                if i < N0-1:
                    if cost_structure == 0:
                        # cost_structure is only time dependent
                        for k in range(K):
                            cb[i, ii, t, k] = fuel * d_II[i,ii] * (1 + costScalingFactor * t)
                    if cost_structure == -1 or cost_structure == 1:
                        for k in range(K):
                            surgeFlag = False;
                            if (t,k) in hurricaneDataSet.nodeTime2Go:
                                if hurricaneDataSet.nodeTime2Go[(t,k)] <= safe_time + 1e-5:
                                    surgeFlag = True;
                            else:
                                surgeFlag = True;
                            if not surgeFlag:
                                cb[i, ii, t, k] = fuel * d_II[i,ii]
                            else:
                                cb[i, ii, t, k] = fuel * d_II[i,ii]*costScalingFactor
                else:
                    if cost_structure == 0:
                        # cost_structure is only time dependent
                        for k in range(K):
                            cb[i, ii, t, k] = fuel * d_KI[ii] * (1 + costScalingFactor * t)
                    if cost_structure == -1 or cost_structure == 1:
                        for k in range(K):
                            surgeFlag = False;
                            if (t,k) in hurricaneDataSet.nodeTime2Go:
                                if hurricaneDataSet.nodeTime2Go[(t,k)] <= safe_time + 1e-5:
                                    surgeFlag = True;
                            else:
                                surgeFlag = True;
                            if not surgeFlag:
                                cb[i, ii, t, k] = fuel * d_KI[ii]
                            else:
                                cb[i, ii, t, k] = fuel * d_KI[ii]*costScalingFactor
    # Unit cost of transporting items from MDC/SP i to/between a demand point j
    ca = np.empty((N0, Nj, T, K))
    for i in range(N0):
        for j in range(Nj):
            for t in range(T):
                if i < N0-1:
                    if cost_structure == 0:
                        # cost_structure is only time dependent
                        for k in range(K):
                            ca[i, j, t, k] = fuel * d_JI[j,i] * (1 + costScalingFactor * t)
                    if cost_structure == -1 or cost_structure == 1:
                        for k in range(K):
                            surgeFlag = False;
                            if (t,k) in hurricaneDataSet.nodeTime2Go:
                                if hurricaneDataSet.nodeTime2Go[(t,k)] <= safe_time + 1e-5:
                                    surgeFlag = True;
                            else:
                                surgeFlag = True;
                            if not surgeFlag:
                                ca[i, j, t, k] = fuel * d_JI[j,i]
                            else:
                                ca[i, j, t, k] = fuel * d_JI[j,i]*costScalingFactor
                else:
                    if cost_structure == 0:
                        # cost_structure is only time dependent
                        for k in range(K):
                            ca[i, j, t, k] = fuel * d_KJ[j] * (1 + costScalingFactor * t)
                    if cost_structure == -1 or cost_structure == 1:
                        for k in range(K):
                            surgeFlag = False;
                            if (t,k) in hurricaneDataSet.nodeTime2Go:
                                if hurricaneDataSet.nodeTime2Go[(t,k)] <= safe_time + 1e-5:
                                    surgeFlag = True;
                            else:
                                surgeFlag = True;
                            if not surgeFlag:
                                ca[i, j, t, k] = fuel * d_KJ[j]
                            else:
                                ca[i, j, t, k] = fuel * d_KJ[j]*costScalingFactor

    cp = np.empty((T, K))
    ch = np.empty((Ni, T))
    for t in range(1, T + 1):
        if cost_structure == 0:
            # cost_structure is only time dependent
            for k in range(K):
                cp[t - 1, k] = base * (1 + costScalingFactor * (t - 1))
        if cost_structure == -1 or cost_structure == 1:
            for k in range(K):
                surgeFlag = False;
                if (t-1,k) in hurricaneDataSet.nodeTime2Go:
                    if hurricaneDataSet.nodeTime2Go[(t-1,k)] <= safe_time + 1e-5:
                        surgeFlag = True;
                else:
                    surgeFlag = True;
                if not surgeFlag:
                    cp[t - 1, k] = base
                else:
                    cp[t - 1, k] = base*costScalingFactor
        ch[:, t - 1] = np.full(Ni, invCostRatio * base)

    p = penCostRatio * base
    q = salvageCostRatio * base
    x_0 = np.zeros(Ni)

    # Demand data
    SCEN = []

    # read Hurricane Florence point forecast
    pf = pd.read_csv(netFolderPath+'Florence_forecast.csv');

    # read some data for the necessary calculation below
    Na = 6; # WARNING: hardcode here!
    aux = pd.read_excel(netFolderPath+'hurricane-position.xlsx');
    study_line = LineString([(aux['line-1-x'][0],aux['line-1-y'][0]),(aux['line-2-x'][0],aux['line-2-y'][0])]);
    for k in range(1, K + 1):
        scen = np.zeros(Nj)
        if states[k-1][2] == T:
            a = states[k - 1][0] # intensity
            l = states[k - 1][1] # forecast error (along the coastline w.r.t. the point forecast)
            hurr_pos = [pf['Longitude'][T-1] + aux['x_rotate'][0]*l/aux['x_miles'][0], pf['Latitude'][T-1] + aux['y_rotate'][0]*l/aux['y_miles'][0]];
            proj_point = study_line.interpolate(study_line.project(Point(hurr_pos)));

            for j in range(1, Nj + 1):
                dLandfall = np.sqrt(pow((df['longitude'][Ni+j]-proj_point.x)*aux['x_miles'][0],2)+pow((df['latitude'][Ni+j]-proj_point.y)*aux['y_miles'][0],2))
                if dLandfall <= cMax:
                    scen[j - 1] = df['Demand'][Ni+j] * (1 - (dLandfall / cMax)) * pow(a-1,2) / pow(Na-1,2)
                else:
                    scen[j - 1] = 0
        SCEN.append(scen)
    x_cap = np.zeros(Ni);
    for i in range(Ni):
        if df['Capacity'].isnull()[1+i]:
            x_cap[i] = 1e8; # hardcode as infinity
        else:
            x_cap[i] = df['Capacity'][1+i];

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

    print("# of SPs = ", self.Ni);
    print("# of DPs = ", self.Nj);  

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
