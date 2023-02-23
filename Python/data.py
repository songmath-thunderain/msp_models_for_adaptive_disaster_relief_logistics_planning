import numpy as np
import pandas as pd
from numpy.linalg import norm

np.random.seed(2022);

#Ni: number of supply points (without MDC).
#Nj: number of demand points.
#N0: number of supply points including MDC.
#x_cap: the capacity of each supply point.
#cb: unit cost of transporting/rerouting relief items from the MDC/SP i to/between different SP i' 
#ch: unit cost for holding an item at SP i.
#h: unit cost for purchasing a relief item.
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

max_iter = 100000;
stall = 200;
cutviol_maxiter = 100000;
nbhrs = 3;
time_limit = nbhrs * 60 ** 2;
eps = 1e-4;

########################################################################################

N0 = Ni + 1;  # number of supply points + the MDC

# the set of possible locations
L = [(0, 100), (100, 200), (200, 300), (300, 400), (400, 500), (500, 600), (600, 700)];  # discussion on the last state

# probability distributions:
P_intesity = pd.read_csv("JuliaData/intensity.csv").values;  # intensity MC
P_location = pd.read_csv("JuliaData/location.csv").values;  # location MC
P_landfall = pd.read_csv("JuliaData/landfall_7.csv").values;  # landfall MC

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

########################################################################################

# locations of the different facilities:
MDC = [350, 450];
x_low = 0;
x_up = 700;
y_low = 0;
y_up = 100;

nodes = pd.read_csv("JuliaData/nodes.csv");  # read the nodes locations

# list for the coordinates of the
# list for the coordinates of the different supply points
SP = [];
for i in range(1, Ni+1):
    SP.append([nodes[i-1, 0], nodes[i-1, 1]+100]);

# list for the coordinates of the different demand points
DP = [];
for j in range(1, Nj+1):
    DP.append([nodes[j-1, 2], nodes[j-1, 3]]);

# fuel cost
fuel = 0.0038;

# unit cost of transporting/rerouting items from MDC/SP i to/between SP i' 
cb = np.empty((N0, Ni, T))
for i in range(1, N0+1):
    for ii in range(1, Ni+1):
        for t in range(1, T+1):
            if i < N0:
                cb[i-1, ii-1, t-1] = fuel*np.linalg.norm(np.array(SP[i-1])-np.array(SP[ii-1]), 2)*(1+factor*(t-1));
            else:
                cb[i-1, ii-1, t-1] = fuel*np.linalg.norm(np.array(MDC)-np.array(SP[ii-1]), 2)*(1+factor*(t-1));

# unit cost of transporting items from MDC/SP i to/between a demand point j
ca = np.empty((N0, Nj, T));
for i in range(1, N0+1):
    for j in range(1, Nj+1):
        for t in range(1, T+1):
            if i < N0:
                ca[i-1, j-1, t-1] = fuel*np.linalg.norm(np.array(SP[i-1])-np.array(DP[j-1]), 2)*(1+factor*(t-1));
            else:
                ca[i-1, j-1, t-1] = fuel*np.linalg.norm(np.array(MDC)-np.array(DP[j-1]), 2)*(1+factor*(t-1));

base = 5;  # base unit cost for logistic costs
h = np.empty(T);
ch = np.empty((Ni, T));
for t in range(1, T+1):
    h[t-1] = base*(1+factor*(t-1));
    ch[:, t-1] = 0.05*base;

p = 50*base;
q = -0.2*base;
D_max = 400;  # the maximum demand that can happen

x_cap = nodes[:Ni, 4]*(Nj/Ni);  # capacity of each SP
x_0 = np.zeros(Ni);  # initial items at different SPs
f_cap = np.full((N0, Ni), np.inf);


# Demand data

SCEN = [];
c_max = 300;  # the largest possible radius of a hurricane
for k in range(K):
    # initialize a scenario matrix (scen) for each DP j and each possible point m in layer2
    scen = np.zeros(Nj); # we will create a scenario list for each state k
    a = S[k][0]; # what is the observed intensity state
    l = S[k][1]; # what is the observed location state
    
    # what are the coordinates of where the hurricane made landfall [xx,yy]
    # note that yy=0 since the hurricane makes landfall in the cost by assumption
    predicted = 0.5 * (L[l][0] + L[l][1]); # this is the predicted x_coordinates for the landfall location: for simplicity, just use the center of the interval	
    xx_coord = max(x_low, min(predicted, x_up)); # we already know x can't be smaller than x_low and bigger than x_up
    landfall = [xx_coord, 0];
    
    # now let's calculate the demand from each DP to the m location 
    for j in range(Nj):
        # how far is the destination DPj j
        c_j = np.linalg.norm(landfall - DP[j]);
        if c_j <= c_max:
            scen[j] = D_max * (1 - (c_j / c_max)) * (a - 1) ** 2 / ((Na - 1) ** 2);
        else:
            scen[j] = 0;
    
    SCEN[k] = scen;
