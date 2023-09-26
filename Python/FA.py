import time
import gurobipy as gp
from gurobipy import GRB
import pandas as pd
import numpy as np
from misc import *
import sys

# Define stage-t problem
def stage_t_state_k_problem(networkDataSet,hurricaneDataSet,t):
    # Create a new model
    m = gp.Model()

    # Data instantiation
    Ni = networkDataSet.Ni;
    N0 = networkDataSet.N0;
    Nj = networkDataSet.Nj;
    T = hurricaneDataSet.T;
    ca = networkDataSet.ca;
    cb = networkDataSet.cb;
    ch = networkDataSet.ch;
    cp = networkDataSet.cp;
    p = networkDataSet.p;
    q = networkDataSet.q;
    x_0 = networkDataSet.x_0;
    x_cap = networkDataSet.x_cap;

    x = {}
    f = {}
    y = {}
    z = {}
    v = {}
    theta = m.addVar(lb = 0)

    for i in range(Ni):
        x[i] = m.addVar(lb=0, ub=x_cap[i])
        v[i] = m.addVar(lb=0)
        for ii in range(Ni):
            f[i, ii] = m.addVar(lb=0)
        for j in range(Nj):
            y[i, j] = m.addVar(lb=0)

    for j in range(Nj):
        z[j] = m.addVar(lb=0)

    # Set objective
    m.setObjective(
        gp.quicksum(cb[i,ii,t] * f[i, ii] for i in range(N0) for ii in range(Ni))
        + gp.quicksum(ch[i,t] * x[i] for i in range(Ni))
        + gp.quicksum(f[N0,i] for i in range(Ni)) * cp[t]
        + gp.quicksum(ca[i,j,t] * y[i, j] for i in range(Ni) for j in range(Nj))
        + gp.quicksum(z[j] for j in range(Nj)) * p
        + gp.quicksum(v[i] for i in range(Ni)) * q
        + theta,
        GRB.MINIMIZE,
    )

    # Define constraints
    FB1Cons = {}  # A dictionary to store flow-balance constraints 1
    FB2Cons = {}  # A dictionary to store flow-balance constraints 2
    dCons = {}  # A dictionary to store all the demand constraints

    for i in range(Ni):
        if t == 0:
            FB1Cons[i] = m.addConstr(
                x[i]
                + gp.quicksum(f[i, j] for j in range(Ni) if j != i)
                - gp.quicksum(f[j, i] for j in range(N0) if j != i)
                + gp.quicksum(y[i, j] for j in range(Nj))
                + v[i]
                == x_0[i]
            )
            FB2Cons[i] = m.addConstr(
                gp.quicksum(f[i, j] for j in range(Ni) if j != i) <= x_0[i]
            )
        else:
            FB1Cons[i] = m.addConstr(
                x[i]
                + gp.quicksum(f[i, j] for j in range(Ni) if j != i)
                - gp.quicksum(f[j, i] for j in range(N0) if j != i)
                + gp.quicksum(y[i, j] for j in range(Nj))
                + v[i]
                == 0
            )
            FB2Cons[i] = m.addConstr(
                gp.quicksum(f[i, j] for j in range(Ni) if j != i) <= 0
            )

    for j in range(Nj):
        dCons[j] = m.addConstr(z[j] + gp.quicksum(y[i, j] for i in range(Ni)) >= 0)

    return m, x, f, y, z, v, theta, dCons, FB1Cons, FB2Cons

# Define the model
def define_models(networkDataSet,hurricaneDataSet,inputParams):
     # Data instantiation
    Ni = networkDataSet.Ni;
    N0 = networkDataSet.N0;
    Nj = networkDataSet.Nj;
    T = hurricaneDataSet.T;
    m = {}
    x = {}
    f = {}
    y = {}
    z = {}
    v = {}
    theta = {}
    dCons = {}
    FB1Cons = {}
    FB2Cons = {}
    for t in range(T):
        for k in hurricaneDataSet.nodeLists[t]:
            ind = k
            (
                m[t, ind],
                x[t, ind],
                f[t, ind],
                y[t, ind],
                z[t, ind],
                v[t, ind],
                theta[t, ind],
                dCons[t, ind],
                FB1Cons[t, ind],
                FB2Cons[t, ind]
            ) = stage_t_state_k_problem(networkDataSet,hurricaneDataSet,t)
            if ind in hurricaneDataSet.absorbing_states:
                theta[t, ind].setAttr(GRB.Attr.UB, 0);
                if inputParams.absorbing_option == 0:
                    for i in range(Ni):
                        for j in range(Ni):
                            if j != i:
                                f[t, ind][i,j].setAttr(GRB.Attr.UB, 0);
                        for j in range(N0):
                            if j != i:
                                f[t, ind][j,i].setAttr(GRB.Attr.UB, 0);
    return m, x, f, y, z, v, theta, dCons, FB1Cons, FB2Cons

# Train model: forward pass
def FOSDDP_forward_pass_oneSP_iteration(networkDataSet,hurricaneDataSet,inputParams, m, x, theta):
    k_init = inputParams.k_init;
    Ni = networkDataSet.Ni;
    T = hurricaneDataSet.T;
    absorbing_states = hurricaneDataSet.absorbing_states;
    k_t = k_init-1;
    in_sample = [k_t];
    xval = [[0] * T for _ in range(Ni)];
    thetaval = [0]*T;
    lb = 1e10;
    for t in range(T):
        if t > 0:
            k_t = MC_sample(in_sample[t - 1], hurricaneDataSet)
            in_sample.append(k_t)
            MSP_fa_update_RHS(m, k_t, t, xval)
        m[t, k_t].optimize()
        if m[t, k_t].status != GRB.OPTIMAL:
            print("Error in Forward Pass")
            print(f"Model in stage = {t} and state = {k_t}, in forward pass is {m[t, k_t].status}")
            sys.exit(0)
        else:
            xval[:, t] = [var.x for var in x[t, k_t]]
            thetaval[t] = theta[t, k_t].x
            if t == 0:
                lb = m[t, k_t].objVal
        if k_t in absorbing_states:
            break
    return xval, thetaval, lb, in_sample


# Train model: backward pass
def FOSDDP_backward_pass_oneSP_iteration(networkDataSet,hurricaneDataSet,inputParams,solveParams, m, x, theta, FB1Cons, FB2Cons, xval, thetaval, in_sample):
    T = hurricaneDataSet.T;
    Ni = networkDataSet.Ni;
    Na = hurricaneDataSet.Na;
    Nb = hurricaneDataSet.Nb;
    K = Na*Nb*T;
    nodeLists = hurricaneDataSet.nodeLists;
    absorbing_states = hurricaneDataSet.absorbing_states;
    cutviolFlag = False;
    for t in range(len(in_sample)-1, 0, -1):
        # Solving all stage-t problems
        Q = [0] * K  #list for all the optimal values
        #list for all the dual multiplies of the first and second set of constraints
        pi1 = [[0] * Ni for _ in range(K)]
        pi2 = [[0] * Ni for _ in range(K)]
        sample_n = in_sample[t-1]; # the state observed at time t-1
        for k in range(len(nodeLists[t])):
            MSP_fa_update_RHS(m, nodeLists[t][k], t, xval)
            m[t, nodeLists[t][k]].optimize()
            if m[t, nodeLists[t][k]].status != GRB.OPTIMAL:
                print("Error in Backward Pass")
                print(f"Model in stage = {t} and state = {nodeLists[t][k]}, in backward pass is {status}")
                sys.exit(0)
            else:
                Q[nodeLists[t][k]] = m[t, nodeLists[t][k]].objVal
                for i in range(Ni):
                    pi1[nodeLists[t][k]][i] = FB1Cons[t,nodeLists[t][k]][i].pi
                    pi2[nodeLists[t][k]][i] = FB2Cons[t,nodeLists[t][k]][i].pi

        # Solving all stage-(t-1) problems and generate cuts/valid inequalities
        for n in range(len(nodeLists[t - 1])):
            if nodeLists[t - 1][n] not in absorbing_states:
                Qvalue = 0
                for k in range(len(nodeLists[t])):
                    if P_joint[nodeLists[t - 1][n]][nodeLists[t][k]] > hurricaneDataSet.smallestTransProb:
                        Qvalue += Q[nodeLists[t][k]] * P_joint[nodeLists[t - 1][n]][nodeLists[t][k]]              					# check if cut is violated at the sample path encountered in the forward pass
                
                # check if cut is violated at the sample path encountered in the forward pass
                if nodeLists[t - 1][n] == sample_n and (
                        (Qvalue - thetaval[t - 1]) / max(1e-10, abs(thetaval[t - 1])) > solveParams.cutviol
                        and abs(Qvalue - thetaval[t - 1]) > solveParams.cutviol
                ):
                    cutviolFlag = True;
                cutcoef = [0] * Ni
                cutrhs_xval = 0
                for k in range(len(nodeLists[t])):
                    if P_joint[nodeLists[t - 1][n]][nodeLists[t][k]] > hurricaneDataSet.smallestTransProb:
                        for i in range(Ni):
                            tempval = (pi1[nodeLists[t][k]][i]
                                        + pi2[nodeLists[t][k]][i]) \
                                        * P_joint[nodeLists[t - 1][n]][nodeLists[t][k]]
                            cutcoef[i] += tempval
                            cutrhs_xval += tempval * xval[i][t - 1]
                m[t - 1, nodeLists[t - 1][n]].addConstr(
                    theta[t - 1, nodeLists[t - 1][n]] - gp.quicksum(
                        cutcoef[i] * x[t - 1, nodeLists[t - 1][n]][i] for i in range(Ni)
                    ) >= Qvalue - cutrhs_xval
                )
    return cutviolFlag

# Train model
def train_models_offline(networkDataSet,hurricaneDataSet,inputParams,solveParams, m, x, theta, FB1Cons, FB2Cons):
    x_0 = networkDataSet.x_0;
    T = hurricaneDataSet.T;
    Ni = networkDataSet.Ni;
    
    # Set the RHS of the first_stage problem
    for i in range(0, Ni):
        FB1Cons[0, i].setAttr(GRB.Attr.RHS, x_0[i]);
        FB2Cons[0, i].setAttr(GRB.Attr.RHS, x_0[i]);

    # Initialize stuff
    train_time = 0
    relative_gap = 1e10
    lb = 0
    LB = []
    xval = [[0] * T for _ in range(Ni)]
    thetaval = [0] * T
    iter = 0
    cutviol_iter = 0
    start = time.time()
    while True:
        iter += 1
        # Forward pass
        xval, thetaval, lb, in_sample = FOSDDP_forward_pass_oneSP_iteration(networkDataSet,hurricaneDataSet,inputParams, m, x, theta)
        LB.append(lb)
        # Termination check
        flag, Elapsed = termination_check(iter, relative_gap, LB, start, cutviol_iter, solveParams)
        if flag != 0:
            train_time = Elapsed
            break
        # Backward pass (if not terminated)
        cutviolFlag = FOSDDP_backward_pass_oneSP_iteration(networkDataSet,hurricaneDataSet,inputParams,solveParams, m, x, theta, FB1Cons, FB2Cons, xval, thetaval, in_sample)
        if cutviolFlag:
            cutviol_iter = 0
        else:
            cutviol_iter += 1
    return LB, train_time, iter

'''
# Evaluate model
def FOSDDP_eval_offline():
    start = time.time()
    osfname = "./data/OOS" + str(k_init) + ".csv"
    OS_paths = np.genfromtxt(osfname, delimiter=",", dtype=int)  # Read the out-of-sample file
    objs_fa = np.zeros((nbOS, T))
    xval_fa = np.empty((nbOS, T), dtype=object)
    fval_fa = np.empty((nbOS, T), dtype=object)
    yval_fa = np.empty((nbOS, T), dtype=object)
    zval_fa = np.empty((nbOS, T), dtype=object)
    vval_fa = np.empty((nbOS, T), dtype=object)
    for s in range(nbOS):
        xval = np.zeros((Ni, T))
        for t in range(1, T + 1):
            k_t = OS_paths[s, t - 1]
            if t > 1:
                MSP_fa_update_RHS(k_t, t, xval)
            m_fa_t_k = m_fa[t, k_t]
            m_fa_t_k.optimize()
            status = m_fa_t_k.status
            if status != GRB.OPTIMAL:
                print(" in evaluation")
                print(f"Model in stage = {t} and state = {k_t}, in forward pass is {status}")
                exit(0)
            else:
                xval_fa[s, t - 1] = [var.x for var in x_fa[t, k_t]]
                xval[:, t - 1] = xval_fa[s, t - 1]
                fval_fa[s, t - 1] = [var.x for var in f_fa[t, k_t]]
                yval_fa[s, t - 1] = [var.x for var in y_fa[t, k_t]]
                zval_fa[s, t - 1] = [var.x for var in z_fa[t, k_t]]
                vval_fa[s, t - 1] = [var.x for var in v_fa[t, k_t]]
                objs_fa[s, t - 1] = m_fa_t_k.objVal - ϴ_fa[t, k_t].x
            if absorbing_option == 0:
                if k_t in absorbing_states:
                    if sum(fval_fa[s, t - 1][N0 - 1][i - 1] for i in range(1, Ni + 1)) > 1e-5:
                        print("Something is wrong! Sum of flow from MDC = ",
                              sum(fval_fa[s, t - 1][N0 - 1][i - 1] for i in range(1, Ni + 1)))
                        exit(0)
            if k_t in absorbing_states:
                break
    fa_bar = np.mean(np.sum(objs_fa, axis=1))
    fa_std = np.std(np.sum(objs_fa, axis=1))
    fa_low = fa_bar - 1.96 * fa_std / np.sqrt(nbOS)
    fa_high = fa_bar + 1.96 * fa_std / np.sqrt(nbOS)
    print("FA...")
    print(f"μ ± 1.96*σ/√NS = {fa_bar} ± {fa_low,fa_high}")
    elapsed = time.time() - start
    vals = [xval_fa, fval_fa, yval_fa, zval_fa, vval_fa]
    return objs_fa, fa_bar, fa_low, fa_high, elapsed

# Update RHS of flow-balance and demand constraint
def MSP_fa_update_RHS(k_t, t, xval):
    for i in range(1, Ni + 1):
        FB1Cons_fa[t - 1][i - 1].RHS = xval[i - 1][t - 2]
        FB2Cons_fa[t - 1][i - 1].RHS = xval[i - 1][t - 2]
    for j in range(1, Nj + 1):
        if S[k_t - 1][2] == Nc and S[k_t - 1][0] != 1:
            dCons_fa[t - 1][j - 1].RHS = SCEN[k_t - 1][j - 1]
        else:
            dCons_fa[t - 1][j - 1].RHS = 0


'''