import numpy as np;
import pandas as pd;
from dataClass import solveParams, hurricaneData, networkData;

def solveParamsInput():
    max_iter = 100000;
    stall = 200;
    cutviol_maxiter = 100000;
    nbhrs = 3;
    time_limit = nbhrs * 60 ** 2;
    eps = 1e-4;
    solveParamSet = solveParams(max_iter,stall,cutviol_maxiter,time_limit,eps);
    return solveParamSet;

def hurricaneInput(intensityFile, locationFile, landfallFile):
    np.random.seed(2023);
    # the set of possible locations
    L = [(0, 100), (100, 200), (200, 300), (300, 400), (400, 500), (500, 600), (600, 700)];  # discussion on the last state

    # probability distributions:
    P_intensity = pd.read_csv(intensityFile).values;  # intensity MC
    P_location = pd.read_csv(locationFile).values;  # location MC
    P_landfall = pd.read_csv(landfallFile).values;  # landfall MC

    Na = P_intesity.shape[0];  # intensity MC number of states
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
                            P_joint[k1-1, k2-1] = P_intesity[k-1, n-1] * P_location[l-1, m-1] * P_landfall[f-1, j-1];
                S[k1-1] = [k, l, f];
                if k == 1 or f == Nc:
                    absorbing_states.append(k1);

    # Create the transition probability from stage 1 to stage T (applying the C-K equation)
    P_terminals = [];
    for t in range(1, T+1):
        P_terminal = np.linalg.matrix_power(P_joint, (T - t));
        P_terminal_c = np.copy(P_terminal);
        # normalize the probabilities
        for k in range(K):
            P_terminal[k, :] /= np.sum(P_terminal_c[k, :]);
        P_terminals.append(P_terminal);

    # normalize the probabilities
    P_temp = np.copy(P_joint);
    for k in range(K):
        P_joint[k, :] /= np.sum(P_temp[k, :]);
        
    hurricaneDataSet = hurricaneData(P_intensity, P_location, P_landfall, T, Tmin, P_joint, absorbing_states, P_terminals
    