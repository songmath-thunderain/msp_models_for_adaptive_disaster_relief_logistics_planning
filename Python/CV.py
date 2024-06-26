import time
import gurobipy as gp
from gurobipy import GRB
import pandas as pd
import numpy as np
import sys

class CV:
    def __init__(self,inputParams,hurricaneData,networkData):
        self.inputParams = inputParams;
        self.hurricaneData = hurricaneData;
        self.networkData = networkData;

    # Define the model
    def deterministic_model(self, sample_path):
        # Define the model
        m = gp.Model("claivoyant");

        # Data instantiation
        Ni = self.networkData.Ni;
        N0 = self.networkData.N0;
        Nj = self.networkData.Nj;
        T = self.hurricaneData.T;
        ca = self.networkData.ca;
        cb = self.networkData.cb;
        cp = self.networkData.cp;
        ch = self.networkData.ch;
        p = self.networkData.p;
        q = self.networkData.q;
        
        x_0 = self.networkData.x_0;
        x_cap = self.networkData.x_cap;
        f_cap = self.networkData.f_cap;
        forbiddenArcs = self.networkData.forbiddenArcs;

        # Define the variables
        x = {};
        for i in range(Ni):
            if x_cap[i] == 1e8:
                for t in range(T):
                    x[i,t] = m.addVar(lb=0);
            else:
                for t in range(T):
                    x[i,t] = m.addVar(lb=0, ub=x_cap[i]);
        f = {};
        for i in range(N0):
            for ii in range(Ni):
                for t in range(T):
                    f[i,ii,t] = m.addVar(lb=0);
        y = {};
        for i in range(Ni):
            for j in range(Nj):
                for t in range(T):
                    y[i,j,t] = m.addVar(lb=0);
        z = {};
        for j in range(Nj):
            for t in range(T):
                z[j,t] = m.addVar(lb=0);
        v = {};
        for i in range(Ni):
            for t in range(T):
                v[i,t] = m.addVar(lb=0);

        # Define the objective
        m.setObjective(
            gp.quicksum(gp.quicksum(gp.quicksum(cb[i, ii, t, sample_path[t]-1] * f[i, ii, t] for ii in range(Ni))
                            for i in range(N0)) for t in range(T))
            + gp.quicksum(gp.quicksum(ch[i, t] * x[i, t] for i in range(Ni)) for t in range(T))
            + gp.quicksum(gp.quicksum(f[N0-1, i, t] for i in range(Ni)) * cp[t, sample_path[t]-1] for t in range(T))
            + gp.quicksum(gp.quicksum(gp.quicksum(ca[i, j, t, sample_path[t]-1] * y[i, j, t] for j in range(Nj))
                                for i in range(Ni)) for t in range(T))
            + gp.quicksum(gp.quicksum(z[j, t] for j in range(Nj)) * p for t in range(T))
            + gp.quicksum(gp.quicksum(v[i, t] for i in range(Ni)) * q for t in range(T)),
            GRB.MINIMIZE
        )

        # Define the constraints

        # flow capacity constraints: total flow per period cannot exceed an upper limit
        for t in range(T):
            m.addConstr(gp.quicksum(gp.quicksum(f[i,ii,t] for i in range(N0)) for ii in range(Ni)) <= f_cap);

        dCons = {}
        for t in range(T):
            for i in range(Ni):
                if t == 0:
                    m.addConstr(
                        x[i, t] + gp.quicksum(f[i, j, t] for j in range(Ni) if j != i)
                        - gp.quicksum(f[j, i, t] for j in range(N0) if j != i)
                        + gp.quicksum(y[i, j, t] for j in range(Nj))
                        + v[i, t] == x_0[i]
                    )
                    m.addConstr(
                        gp.quicksum(f[i, j, t] for j in range(Ni) if j != i) <= x_0[i]
                    )
                else:
                    m.addConstr(
                        x[i, t] + gp.quicksum(f[i, j, t] for j in range(Ni) if j != i)
                        - gp.quicksum(f[j, i, t] for j in range(N0) if j != i)
                        + gp.quicksum(y[i, j, t] for j in range(Nj))
                        + v[i, t] == x[i, t-1]
                    )
                    m.addConstr(gp.quicksum(f[i,j,t] for j in range(Ni) if j != i) <= x[i,t-1])
            for j in range(Nj):
                dCons[t,j] = m.addConstr(z[j,t]+gp.quicksum(y[i,j,t] for i in range(Ni)) >= 0)

        for i in range(Ni):
            for j in range(Nj):
                if (i,j) in forbiddenArcs:
                    for t in range(T):
                        y[i,j,t].setAttr(GRB.Attr.UB, 0);

        m.update();
        m.setParam("OutputFlag", 0);

        return m, x, f, y, z, v, dCons;
    
    
    def clairvoyant_eval(self,osfname):
        start = time.time()

        absorbing_option = self.inputParams.absorbing_option;
        nbOS = self.inputParams.nbOS;

        T = self.hurricaneData.T;
        Ni = self.networkData.Ni;
        N0 = Ni + 1;
        Nj = self.networkData.Nj;
        SCEN = self.networkData.SCEN;
        S = self.hurricaneData.states;

        OS_paths = pd.read_csv(osfname).values  # Read the out-of-sample file

        objs_cv = np.zeros(nbOS)

        for s in range(nbOS):
            absorbingT = list(OS_paths[s, 0:T]).index(next((x for x in OS_paths[s, 0:T] if (x-1) in self.hurricaneData.absorbing_states), None))

            if S[OS_paths[s, absorbingT]-1][0] == 1:
                continue
            else:
                # Define the clairvoyant model
                m_cv, x_cv, f_cv, y_cv, z_cv, v_cv, dCons_cv = self.deterministic_model(OS_paths[s, 0:T])
                k_t = OS_paths[s, absorbingT]  # State corresponding to the landfall
                
                for j in range(Nj):
                    for t in range(T):
                        if t == absorbingT:
                            dCons_cv[absorbingT,j].setAttr(GRB.Attr.RHS, SCEN[k_t-1][j]) # Set the RHS of demand constraint
                            if absorbing_option == 0:
                                for i in range(N0):
                                    for ii in range(Ni):
                                        f_cv[i,ii,t].setAttr(GRB.Attr.UB, 0);
                        else:
                            dCons_cv[t,j].setAttr(GRB.Attr.RHS, 0);

                m_cv.optimize()  # Solve the model
                if m_cv.status != GRB.status.OPTIMAL:
                    print(" Model in clairvoyant is ", m_cv.status)
                    sys.exit(0)
                else:
                    #xval_cv[s] = [var.x for var in x_cv]
                    #fval_cv[s] = [var.x for var in f_cv]
                    #yval_cv[s] = [var.x for var in y_cv]
                    #zval_cv[s] = [var.x for var in z_cv]
                    #vval_cv[s] = [var.x for var in v_cv]
                    objs_cv[s] = m_cv.ObjVal;

        cv_bar = np.mean(objs_cv)
        cv_std = np.std(objs_cv)
        cv_low = cv_bar - 1.96 * cv_std / np.sqrt(nbOS)
        cv_high = cv_bar + 1.96 * cv_std / np.sqrt(nbOS)
        CI = cv_bar - cv_low;
        print("Clairvoyant....")
        print(f"μ ± 1.96*σ/√NS = {cv_bar} ± {CI}")

        elapsed_cv = time.time() - start
        return [cv_bar, CI, elapsed_cv]

