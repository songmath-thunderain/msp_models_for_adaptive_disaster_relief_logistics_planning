import time
import gurobipy as gp
from gurobipy import GRB
import pandas as pd
import numpy as np
from misc import *
import sys

# Define the terminal-stage problem: only used when absorbing_option = 1, i.e., MDC/SP operation is allowed to occur
def terminal_model(networkDataSet,hurricaneDataSet, t_roll, x_init):
    # Define the model
    m = gp.Model()
    # Data instantiation
    Ni = networkDataSet.Ni;
    N0 = networkDataSet.N0;
    Nj = networkDataSet.Nj;
    ca = networkDataSet.ca;
    cb = networkDataSet.cb;
    ch = networkDataSet.ch;
    cp = networkDataSet.cp;
    p = networkDataSet.p;
    q = networkDataSet.q;
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
        for j in range(Nj):
            y[i, j] = m.addVar(lb=0)
    
    for i in range(N0):
        for ii in range(Ni):
            f[i, ii] = m.addVar(lb=0)

    for j in range(Nj):
        z[j] = m.addVar(lb=0)
    
    # Define the objective
    m.setObjective(
        gp.quicksum(cb[i,ii,t_roll] * f[i, ii] for ii in range(Ni) for i in range(N0))
        + gp.quicksum(ch[i,t_roll] * x[i] for i in range(Ni))
        + gp.quicksum(f[N0-1, i] for i in range(Ni)) * cp[t_roll]
        + gp.quicksum(ca[i,j,t_roll] * y[i, j] for i in range(Ni) for j in range(Nj))
        + gp.quicksum(z[j] for j in range(Nj)) * p
        + gp.quicksum(v[i] for i in range(Ni)) * q,
        GRB.MINIMIZE)
    
    # Define the constraints
    dCons = {}
    for i in range(Ni):
        m.addConstr(
            x[i] + gp.quicksum(f[i, j] for j in range(Ni) if j != i)
            - gp.quicksum(f[j, i] for j in range(N0) if j != i)
            + gp.quicksum(y[i, j] for j in range(Nj))
            + v[i] == x_init[i]
        )
        m.addConstr(gp.quicksum(f[i, j] for j in range(Ni) if j != i) <= x_init[i])
    
    for j in range(Nj):
        dCons[j] = m.addConstr(z[j] + gp.quicksum(y[i, j] for i in range(Ni)) >= 0)

    m.update();
    m.setParam("OutputFlag", 0);
    return m, x, f, y, z, v, dCons

# Define first-stage master problem
def RH_2SSP_first_stage(networkDataSet, hurricaneDataSet, inputParams, t_roll, k_t, x_init):
    # Note that the static 2SSP corresponds to the case when t_roll = 0
    # Data instantiation
    Ni = networkDataSet.Ni;
    N0 = networkDataSet.N0;
    Nj = networkDataSet.Nj;
    ca = networkDataSet.ca;
    cb = networkDataSet.cb;
    ch = networkDataSet.ch;
    cp = networkDataSet.cp;
    p = networkDataSet.p;
    q = networkDataSet.q;
    x_cap = networkDataSet.x_cap;
    T = hurricaneDataSet.T;
    absorbing_option = inputParams.absorbing_option;
    nodeScenList = hurricaneDataSet.nodeScenList;
    nodeScenWeights = hurricaneDataSet.nodeScenWeights;

    nbstages1 = T - t_roll
    if absorbing_option == 0:
        # If no MDC/SP operation is allowed in the absorbing state, do not plan for stage T
        nbstages1 = T - t_roll - 1
    
    # Define the model
    m = gp.Model()
    
    nbScens = len(nodeScenList[t_roll][k_t])
    pbScens = nodeScenWeights[t_roll][k_t]
    
    # Define the variables
    x = {}
    f = {}
    theta = {}

    for t in range(nbstages1):
        for i in range(Ni):
            x[i,t] = m.addVar(lb=0, ub=x_cap[i])
    
        for i in range(N0):
            for ii in range(Ni):
                f[i, ii, t] = m.addVar(lb=0)

    for n in range(nbScens):
        theta[n] = m.addVars(lb=-1e8)
    
    # Define the objective
    m.setObjective(
        gp.quicksum(
            gp.quicksum(
                gp.quicksum(cb[i,ii,t+t_roll] * f[i, ii, t] for ii in range(Ni) for i in range(N0))
                for t in range(nbstages1))
            for i in range(Ni))
        + gp.quicksum(
            gp.quicksum(ch[i,t_roll+t] * x[i, t] for i in range(Ni))
            for t in range(nbstages1))
        + gp.quicksum(
            gp.quicksum(cp[t_roll+t] * f[N0, i, t] for i in range(Ni))
            for t in range(nbstages1))
        + gp.quicksum(theta[n] * pbScens[n] for n in range(nbScens)),
        GRB.MINIMIZE)
    
    # Define the constraints
    for t in range(nbstages1):
        for i in range(Ni):
            if t == 0:
                m.addConstr(
                    x[i, t] + gp.quicksum(f[i, j, t] for j in range(Ni) if j != i)
                    - gp.quicksum(f[j, i, t] for j in range(N0) if j != i) == x_init[i]
                )
                m.addConstr(gp.quicksum(f[i, j, t] for j in range(Ni) if j != i) <= x_init[i])
            else:
                m.addConstr(
                    x[i, t-1] - x[i, t] - gp.quicksum(f[i, j, t] for j in range(Ni) if j != i)
                    + gp.quicksum(f[j, i, t] for j in range(N0) if j != i) == 0
                )
                m.addConstr(gp.quicksum(f[i, j, t] for j in range(Ni) if j != i) <= x[i, t-1])
    
    return m, x, f, theta

# Define second-stage scenario subproblem
def RH_2SSP_second_stage(networkDataSet):
    # Data instantiation
    Ni = networkDataSet.Ni;
    N0 = networkDataSet.N0;
    Nj = networkDataSet.Nj;
    ca = networkDataSet.ca;
    cb = networkDataSet.cb;
    ch = networkDataSet.ch;
    cp = networkDataSet.cp;
    p = networkDataSet.p;
    q = networkDataSet.q;
    # Define the model
    m = gp.Model()
    
    # Define the variables
    y = {}
    z = {}
    v = {}
    for i in range(Ni):
        for j in range(Nj):
            y[i,j] = m.addVar(lb=0)

    for j in range(Nj):
        z[j] = m.addVar(lb = 0)

    for i in range(Ni):
        v[i] = m.addVar(lb = 0)
    reimbursement = m.addVar(ub=0, lb = -GRB.INFINITY)
    
    # Define the objective
    m.setObjective(
        gp.quicksum(ca[i,j,0] * y[i, j] for i in range(Ni) for j in range(Nj))
        + gp.quicksum(z[j] for j in range(Nj)) * p
        + gp.quicksum(v[i] for i in range(Ni)) * q
        + reimbursement,
        GRB.MINIMIZE)
    
    # Define the constraints
    xCons = {}
    dCons = {}
    
    for i in range(Ni):
        xCons[i] = m.addConstr(gp.quicksum(y[i, j] for j in range(Nj)) + v[i] == 0)
    
    for j in range(Nj):
        dCons[j] = m.addConstr(z[j] + gp.quicksum(y[i, j] for i in range(Ni)) == 0)
    
    rCons = m.addConstr(reimbursement == 0)
    
    return m, y, xCons, dCons, rCons


# Define the two-stage SP models: master problem and subproblem
def RH_2SSP_define_models(networkDataSet, hurricaneDataSet, inputParams, t_roll, k_t, x_init):
    # Define first stage (master problem) model
    master, x, f, theta = RH_2SSP_first_stage(networkDataSet, hurricaneDataSet, inputParams, t_roll, k_t, x_init)
    
    # Define second stage (subproblem) optimality model
    subproblem, y2, xCons, dCons, rCons = RH_2SSP_second_stage(networkDataSet)
    
    return master, x, f, theta, subproblem, y2, xCons, dCons, rCons


# Solves the two-stage SP model
def RH_2SSP_solve_roll(networkDataSet, hurricaneDataSet, inputParams, solveParams, k_t, t_roll, master, subproblem, x, f, theta, y, xCons, dCons, rCons):
    Ni = networkDataSet.Ni;
    N0 = networkDataSet.N0;
    T = hurricaneDataSet.T;
    absorbing_option = inputParams.absorbing_option;
    nodeScenList = hurricaneDataSet.nodeScenList;
    nodeScenWeights = hurricaneDataSet.nodeScenWeights;

    nbstages1 = T - t_roll
    if absorbing_option == 0:
        # If no MDC/SP operation is allowed in the absorbing state, do not plan for stage T since we know for sure that all states are absorbing
        nbstages1 = T - t_roll - 1
    nbScens = len(nodeScenList[t_roll][k_t])
    pbScens = nodeScenWeights[t_roll][k_t]
    LB = -1e10
    UB = 1e10
    thetaval = [0] * nbScens
    xval = [[0] * nbstages1 for _ in range(Ni)]
    fval = [[[0] * nbstages1 for _ in range(Ni)] for _ in range(N0)]
    solveIter = 0
    flag = 1
    
    while (UB - LB) * 1.0 / max(1e-10, abs(LB)) > solveParams.cutviol and abs(UB - LB) > solveParams.cutviol and solveIter < 100:  # WARNING: Temporary fix here!
        # Solve first stage
        flag = 0
        solveIter += 1
        LB, xval, fval, thetaval = solve_first_stage(LB, xval, fval, thetaval, master, x, f, theta, nbstages1, nbScens)
        firstCost = LB - sum(thetaval[n] * pbScens[n] for n in range(nbScens))
        
        # Solve second stage
        flag, Qbar = solve_second_stage(k_t, t_roll, xval, fval, thetaval, master, subproblem, x, f, theta, y, xCons, dCons, rCons)
        if flag != -1:
            UB = min(firstCost + Qbar, UB)
    
    if solveIter == 100:
        print("# iterations is maxed out!")
        print("LB =", LB)
        print("UB =", UB)
    
    return LB, UB, xval, fval, thetaval

# Solves the first-stage problem
def solve_first_stage(networkDataSet, master, x, f, theta, nbstages1, nbScens):
    Ni = networkDataSet.Ni;
    N0 = networkDataSet.N0;
    master.optimize()
    status_master = master.status
    if status_master != GRB.OPTIMAL:
        print("Master problem status is:", status_master, "=== Oops! :/")
        sys.exit(0)
    else:
        LB = master.objVal
        xval = [[x[i, t].x for t in range(nbstages1)] for i in range(Ni)]
        fval = [[[f[i, j, t].x for t in range(nbstages1)] for j in range(Ni)] for i in range(N0)]
        thetaval = [theta[n].x for n in range(nbScens)]
        return LB, xval, fval, thetaval

# Solves the second-stage problem
def solve_second_stage(networkDataSet, hurricaneDataSet, inputParams, solveParams, k_t, t_roll, xval, fval, thetaval, master, subproblem, x, f, theta, y, xCons, dCons, rCons):
    Ni = networkDataSet.Ni;
    N0 = networkDataSet.N0;
    T = hurricaneDataSet.T;
    cb = networkDataSet.cb;
    ch = networkDataSet.ch;
    cp = networkDataSet.cp;
    absorbing_option = inputParams.absorbing_option;
    nodeScenList = hurricaneDataSet.nodeScenList;
    nodeScenWeights = hurricaneDataSet.nodeScenWeights;
    
    nbScens = len(nodeScenList[t_roll][k_t])
    pbScens = nodeScenWeights[t_roll][k_t]
    flag = 0
    nbstages1 = T - t_roll + 1
    
    if absorbing_option == 0:
        nbstages1 = T - t_roll
    
    Q = [0] * nbScens
    pi1 = [[] for _ in range(nbScens)]
    pi2 = [[] for _ in range(nbScens)]
    pi3 = [0] * nbScens
    Qbar = 0
    
    for n in range(nbScens):
        absorbingT = nodeScenList[t_roll][k_t][n][0]
        RH_2SSP_update_RHS(absorbingT, k_t, subproblem, xCons, dCons, rCons, xval, fval, y, t_roll)
        Q[n], pi1[n], pi2[n], pi3[n], flag = solve_scen_subproblem(subproblem, xCons, dCons, rCons)
        
        if flag == -1:
            print("Subproblem status is infeasible?!")
            exit(0)
    
    Qbar = sum(Q[n] * pbScens[n] for n in range(nbScens))
    
    # Cut generation: multi-cut version
    for n in range(nbScens):
        if (Q[n] - thetaval[n]) / max(1e-10, abs(Q[n])) > solveParams.cutviol and Q[n] - thetaval[n] > solveParams.cutviol:
            tt = nodeScenList[t_roll][k_t][n][0]
            
            if absorbing_option == 0:
                master.addConstr(
                    theta[n] - sum(pi1[n][i] * x[i, tt - t_roll - 1] for i in range(Ni))
                    - pi3[n] * (-sum(
                        sum(
                            cb[i,ii,t_roll + t] * f[i,ii,t]
                            for ii in range(Ni)
                        ) for i in range(N0))
                    + sum(ch[i][t_roll + t] * x[i, t] for i in range(Ni))
                    + sum(f[N0,i,t] for i in range(Ni)) * cp[t_roll + t]
                    for t in range(tt - t_roll, nbstages1)
                    ) >= Q[n] - sum(pi1[n][i] * xval[i,tt - t_roll] for i in range(Ni)) - pi3[n] * (-sum(
                        sum(
                            cb[i,ii,t_roll + t] * fval[i,ii,t]
                            for ii in range(Ni)
                        ) for i in range(N0)))
                    + sum(ch[i,t_roll + t] * xval[i,t] for i in range(Ni))
                    + sum(fval[N0-1,i,t] for i in range(Ni)) * cp[t_roll + t]
                )
            else:
                master.addConstr(
                    theta[n] - sum(pi1[n][i] * x[i, tt - t_roll] for i in range(Ni))
                    - pi3[n] * (-sum(
                        sum(
                            cb[i,ii,t_roll + t] * f[i,ii,t]
                            for ii in range(Ni)
                        ) for i in range(N0))
                    + sum(ch[i,t_roll + t] * x[i, t] for i in range(Ni))
                    + sum(f[N0-1,i,t] for i in range(Ni)) * cp[t_roll + t]
                    for t in range(tt + 1 - t_roll, nbstages1)
                    ) >= Q[n] - sum(pi1[n][i] * xval[i,tt - t_roll] for i in range(Ni)) - pi3[n] * (-sum(
                        sum(
                            cb[i,ii,t_roll + t] * fval[i,ii,t]
                            for ii in range(Ni)
                        ) for i in range(N0)))
                    + sum(ch[i,t_roll + t] * xval[i,t] for i in range(Ni))
                    + sum(fval[N0-1,i,t] for i in range(Ni)) * cp[t_roll + t]
                )
            
            flag = 1
    
    return flag, Qbar

# Solves scenario subproblem of the second stage
def solve_scen_subproblem(networkDataSet, subproblem, xCons, dCons, rCons):
    Ni = networkDataSet.Ni;
    Nj = networkDataSet.Nj;
    subproblem.optimize()
    status_subproblem = subproblem.status
    
    Qtemp = 0
    pi1temp = [0] * Ni
    pi2temp = [0] * Nj
    pi3temp = 0
    
    if status_subproblem != GRB.OPTIMAL:
        if status_subproblem == GRB.INFEASIBLE:
            flag = -1
    else:
        Qtemp = subproblem.objVal
        pi1temp = [xCons[i].pi for i in range(Ni)]
        pi2temp = [dCons[j].pi for j in range(Nj)]
        pi3temp = rCons.pi
    
    return Qtemp, pi1temp, pi2temp, pi3temp, flag

# Updates the RHS of the 2nd-stage constraints and objective coefficients
def RH_2SSP_update_RHS(networkDataSet, hurricaneDataSet, inputParams, absorbingT, k_t, subproblem, xCons, dCons, rCons, xval, fval, y, t_roll):
    Ni = networkDataSet.Ni;
    N0 = networkDataSet.N0;
    SCEN = networkDataSet.SCEN;
    T = hurricaneDataSet.T;
    cb = networkDataSet.cb;
    ch = networkDataSet.ch;
    cp = networkDataSet.cp;
    absorbing_option = inputParams.absorbing_option;
    nbstages1 = T - t_roll + 1
    
    if absorbing_option == 0:
        nbstages1 = T - t_roll
    
    for i in range(Ni):
        if absorbing_option == 0:
            xCons[i].setAttr(GRB.Attr.RHS, xval[i][absorbingT - t_roll - 1])
        else:
            xCons[i].setAttr(GRB.Attr.RHS, xval[i][absorbingT - t_roll])

    for j in range(Nj):
        if S[k_t][0] != 1:
            dCons[j].setAttr(GRB.Attr.RHS, SCEN[k_t][j]);
        else:
            dCons[j].setAttr(GRB.Attr.RHS, 0);
    
    if absorbingT == T:
        rCons.setAttr(GRB.Attr.RHS, 0);
    else:
        updatedRHS = 0
        
        if absorbing_option == 0:
            updatedRHS = -sum(
                sum(
                    cb[i,ii,t_roll + t] * fval[i,ii,t]
                    for ii in range(Ni)
                ) for i in range(N0)
            ) + sum(ch[i,t_roll + t] * xval[i,t] for i in range(Ni)) + sum(fval[N0-1,i,t] for i in range(Ni)) * cp[t_roll + t]
        else:
            updatedRHS = -sum(
                sum(
                    cb[i,ii,t_roll + t] * fval[i,ii,t]
                    for ii in range(Ni)
                ) for i in range(N0)
            ) + sum(ch[i,t_roll + t] * xval[i,t] for i in range(Ni)) + sum(fval[N0-1,i,t] for i in range(Ni)) * cp[t_roll + t]
        
        rCons.setAttr(GRB.Attr.RHS, updatedRHS);
    
    return subproblem

# Main function to solve the rolling horizon problem
def RH_2SSP_main():
    k = 1
    LB = 0
    UB = 1e10
    master, x, f, θ, subproblem, y2, xCons, dCons, rCons = RH_2SSP_define_models(t_roll, k, x_init)
    
    while k <= K:
        print("******")
        print("Rolling horizon iteration k =", k)
        print("******")
        
        # Solve RH_2SSP
        LB, UB, xval, fval, θval = RH_2SSP_solve_roll(k, t_roll, master, subproblem, x, f, θ, y2, xCons, dCons, rCons)
        
        # Update t_roll
        t_roll += 1
        if absorbing_option == 1:
            if t_roll > T:
                t_roll = 1
                k += 1
        else:
            if t_roll > T - 1:
                t_roll = 1
                k += 1
    
    return LB, UB, xval, fval, θval

