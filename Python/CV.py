import time
from gurobipy import *
import pandas as pd
import numpy as np

# Define the model
def deterministic_model():
    # Define the model
    m = Model(optimizer_with_attributes(() -> Gurobi.Optimizer(GRB_ENV), "OutputFlag" => 0))

    # Define the variables
    x = m.addVars(Ni, T, lb=0, ub=x_cap, name="x")
    f = m.addVars(N0, Ni, T, lb=0, ub=f_cap, name="f")
    y = m.addVars(Ni, Nj, T, lb=0, name="y")
    z = m.addVars(Nj, T, lb=0, name="z")
    v = m.addVars(Ni, T, lb=0, name="v")

    # Define the objective
    m.setObjective(
        quicksum(quicksum(quicksum(cb[i, ii, t] * f[i, ii, t] for ii in range(1, Ni+1))
                          for i in range(1, N0+1)) for t in range(1, T+1))
        + quicksum(quicksum(ch[i, t] * x[i, t] for i in range(1, Ni+1)) for t in range(1, T+1))
        + quicksum(quicksum(f[N0, i, t] for i in range(1, Ni+1)) * h[t] for t in range(1, T+1))
        + quicksum(quicksum(quicksum(ca[i, j, t] * y[i, j, t] for j in range(1, Nj+1))
                            for i in range(1, Ni+1)) for t in range(1, T+1))
        + quicksum(quicksum(z[j, t] for j in range(1, Nj+1)) * p for t in range(1, T+1))
        + quicksum(quicksum(v[i, t] for i in range(1, Ni+1)) * q for t in range(1, T+1)),
        GRB.MINIMIZE
    )

    # Define the constraints
    dCons = {}
    for t in range(1, T+1):
        for i in range(1, Ni+1):
            if t == 1:
                m.addConstr(
                    x[i, t] + quicksum(f[i, j, t] for j in range(1, Ni+1) if j != i)
                    - quicksum(f[j, i, t] for j in range(1, N0+1) if j != i)
                    + quicksum(y[i, j, t] for j in range(1, Nj+1))
                    + v[i, t] == x_0[i]
                )
                m.addConstr(
                    quicksum(f[i, j, t] for j in range(1, Ni+1) if j != i) <= x_0[i]
                )
            else:
                m.addConstr(
                    x[i, t] + quicksum(f[i, j, t] for j in range(1, Ni+1) if j != i)
                    - quicksum(f[j, i, t] for j in range(1, N0+1) if j != i)
                    + quicksum(y[i, j, t] for j in range(1, Nj+1))
                    + v[i, t] == x[i, t-1]
                )
                m.addConstr(sum(f[i,j,t] for j in range(1, Ni+1) if j != i) <= x[i,t-1])
        for j in range(1, Nj+1):
            dCons[t,j] = m.addConstr(z[j,t]+sum(y[i,j,t] for i in range(1, Ni+1)) >= 0)

    return m, x, f, y, z, v, dCons;
    
    
def clairvoyant_eval():
    start = time.time()
    osfname = "JuliaData/OOS" + str(k_init) + ".csv"
    OS_paths = pd.read_csv(osfname, header=None).values

    objs_cv = np.zeros(nbOS)
    xval_cv, fval_cv, yval_cv, zval_cv, vval_cv = [], [], [], [], []

    for s in range(nbOS):
        # find the period when the hurricane made landfall && intensity != 1
        tau = np.where((S[OS_paths[s, :T]-1, 2] == Nc-1) & (OS_paths[s, :T] - 1 != absorbing_states))[0]
        if tau.size == 0:
            continue
        else:
            tau = tau[0]
            k_t = OS_paths[s, tau] # state corresponding to the landfall
            for j in range(Nj):
                for t in range(T):
                    if t == tau:
                        dCons_cv[tau, j].setAttr(GRB.Attr.RHS, SCEN[k_t][j])
                    else:
                        dCons_cv[t, j].setAttr(GRB.Attr.RHS, 0)

            m_cv.optimize() # solve the model
            status = m_cv.status # check the status
            if status != GRB.OPTIMAL:
                print("Model in clairvoyant is ", status)
                exit(0)
            else:
                xval_cv.append(m_cv.getAttr('x', x_cv))
                fval_cv.append(m_cv.getAttr('x', f_cv))
                yval_cv.append(m_cv.getAttr('x', y_cv))
                zval_cv.append(m_cv.getAttr('x', z_cv))
                vval_cv.append(m_cv.getAttr('x', v_cv))
                objs_cv[s] = m_cv.objVal

    cv_bar = np.mean(objs_cv)
    cv_std = np.std(objs_cv)
    cv_low = cv_bar - 1.96 * cv_std / np.sqrt(nbOS)
    cv_high = cv_bar + 1.96 * cv_std / np.sqrt(nbOS)
    print("Clairvoyant....")
    print("μ ± 1.96*σ/√NS = {}±{}".format(cv_bar, [cv_low, cv_high]))
    elapsed = time.time() - start
    vals = [xval_cv, fval_cv, yval_cv, zval_cv, vval_cv]

    return objs_cv, cv_bar, cv_low, cv_high, elapsed