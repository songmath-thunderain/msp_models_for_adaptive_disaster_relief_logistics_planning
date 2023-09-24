import time
from gurobipy import *
import pandas as pd
import numpy as np

# Define the model
def deterministic_model(networkDataSet,hurricaneDataSet):
    # Define the model
    m = Model(optimizer_with_attributes(() -> Gurobi.Optimizer(GRB_ENV), "OutputFlag" => 0))

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

    # Define the variables
    x = m.addVars(Ni, T, lb=0, ub=x_cap, name="x")
    f = m.addVars(N0, Ni, T, lb=0, name="f")
    y = m.addVars(Ni, Nj, T, lb=0, name="y")
    z = m.addVars(Nj, T, lb=0, name="z")
    v = m.addVars(Ni, T, lb=0, name="v")

    # Define the objective
    m.setObjective(
        quicksum(quicksum(quicksum(cb[i, ii, t] * f[i, ii, t] for ii in range(1, Ni+1))
                          for i in range(1, N0+1)) for t in range(1, T+1))
        + quicksum(quicksum(ch[i, t] * x[i, t] for i in range(1, Ni+1)) for t in range(1, T+1))
        + quicksum(quicksum(f[N0, i, t] for i in range(1, Ni+1)) * cp[t] for t in range(1, T+1))
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
    
    
def clairvoyant_eval(networkDataSet,hurricaneDataSet,inputParams):
    start = time.time()

    dissipate_option = inputParams.dissipate_option;
    absorbing_option = inputParams.absorbing_option;
    k_init = inputParams.k_init;
    nbOS = inputParams.nbOS;

    T = hurricaneDataSet.T;
    Ni = networkDataSet.Ni;
    Nj = networkDataSet.Nj;
    SCEN = networkDataSet.SCEN;
    S = hurricaneDataSet.states;

    osfname = f"data/synthetic/OOS{k_init}.csv"
    OS_paths = pd.read_csv(osfname).values  # Read the out-of-sample file

    objs_cv = np.zeros(nbOS)
    xval_cv = [None] * nbOS
    fval_cv = [None] * nbOS
    yval_cv = [None] * nbOS
    zval_cv = [None] * nbOS
    vval_cv = [None] * nbOS

    for s in range(nbOS):
        absorbingT = -1

        if dissipate_option:
            absorbingT = list(OS_paths[s, 0:T]).index(next((x for x in OS_paths[s, 0:T] if S[x][2] == T or S[x][0] == 1), None))
        else:
            absorbingT = list(OS_paths[s, 0:T]).index(next((x for x in OS_paths[s, 0:T] if S[x][2] == T), None))

        if OS_paths[s, absorbingT] == 1:
            continue
        else:
            # Define the clairvoyant model
            m_cv, x_cv, f_cv, y_cv, z_cv, v_cv, dCons_cv = deterministic_model(networkDataSet,hurricaneDataSet)
            k_t = OS_paths[s, absorbingT]  # State corresponding to the landfall

            for j in range(Nj):
                for t in range(T):
                    if t == absorbingT:
                        set_normalized_rhs(dCons_cv[absorbingT, j], SCEN[k_t][j])  # Set the RHS of demand constraint

                        if absorbing_option == 0:
                            for i in range(N0):
                                for ii in range(Ni):
                                    set_upper_bound(f_cv[i, ii, t], 0)
                    else:
                        set_normalized_rhs(dCons_cv[t, j], 0)  # Set the RHS of demand constraint

            optimize(m_cv)  # Solve the model
            status = termination_status(m_cv)  # Check the status

            if status != MOI.OPTIMAL:
                print(" Model in clairvoyant is ", status)
                exit(0)
            else:
                xval_cv[s] = [value(var) for var in x_cv]
                fval_cv[s] = [value(var) for var in f_cv]
                yval_cv[s] = [value(var) for var in y_cv]
                zval_cv[s] = [value(var) for var in z_cv]
                vval_cv[s] = [value(var) for var in v_cv]
                objs_cv[s] = objective_value(m_cv)

    cv_bar = np.mean(objs_cv)
    cv_std = np.std(objs_cv)
    cv_low = cv_bar - 1.96 * cv_std / np.sqrt(nbOS)
    cv_high = cv_bar + 1.96 * cv_std / np.sqrt(nbOS)

    print("Clairvoyant....")
    print(f"μ ± 1.96*σ/√NS = {cv_bar} ± [{cv_low}, {cv_high}]")

    elapsed_cv = time.time() - start
    vals = [xval_cv, fval_cv, yval_cv, zval_cv, vval_cv]

    fname = "./output/benchmark/CVresults.csv"
    df = pd.read_csv(fname)
    results_cv = df.values
    results_cv[inst, 0] = 0
    results_cv[inst, 1] = cv_bar
    results_cv[inst, 2] = cv_bar - cv_low
    results_cv[inst, 3] = cv_bar + cv_low
    results_cv[inst, 4] = elapsed_cv
    results_cv[inst, 5] = 0

    updf = pd.DataFrame(results_cv, columns=df.columns)
    updf.to_csv(fname, index=False)

    return objs_cv, cv_bar, cv_low, cv_high, elapsed