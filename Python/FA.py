import time
import gurobipy as gp
from gurobipy import GRB
import pandas as pd
import numpy as np
from misc import *
import sys

class FA:
    def __init__(self,inputParams,solveParams,hurricaneData,networkData):
        self.inputParams = inputParams;
        self.solveParams = solveParams;
        self.hurricaneData = hurricaneData;
        self.networkData = networkData;


    # Define stage-t problem
    def stage_t_state_k_problem(self,t,k):
        # Create a new model
        m = gp.Model()

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

        x = {}
        f = {}
        y = {}
        z = {}
        v = {}
        theta = m.addVar(lb = 0)

        for i in range(Ni):
            if x_cap[i] == 1e8:
                x[i] = m.addVar(lb=0);
            else:
                x[i] = m.addVar(lb=0, ub=x_cap[i])
            v[i] = m.addVar(lb=0)
            for j in range(Nj):
                y[i, j] = m.addVar(lb=0)
        
        for i in range(N0):
            for ii in range(Ni):
                f[i, ii] = m.addVar(lb=0)

        for j in range(Nj):
            z[j] = m.addVar(lb=0)

        # Set objective
        m.setObjective(
            gp.quicksum(cb[i,ii,t,k] * f[i, ii] for i in range(N0) for ii in range(Ni))
            + gp.quicksum(ch[i,t] * x[i] for i in range(Ni))
            + gp.quicksum(f[N0-1,i] for i in range(Ni)) * cp[t,k]
            + gp.quicksum(ca[i,j,t,k] * y[i, j] for i in range(Ni) for j in range(Nj))
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

        m.update();
        m.setParam("OutputFlag", 0);

        return m, x, f, y, z, v, theta, dCons, FB1Cons, FB2Cons

    # Define the model
    def define_models(self):
        # Data instantiation
        Ni = self.networkData.Ni;
        N0 = self.networkData.N0;
        Nj = self.networkData.Nj;
        T = self.hurricaneData.T;
        self.m = {}
        self.x = {}
        self.f = {}
        self.y = {}
        self.z = {}
        self.v = {}
        self.theta = {}
        self.dCons = {}
        self.FB1Cons = {}
        self.FB2Cons = {}
        for t in range(T):
            for k in self.hurricaneData.nodeLists[t]:
                ind = k
                (
                    self.m[t, ind],
                    self.x[t, ind],
                    self.f[t, ind],
                    self.y[t, ind],
                    self.z[t, ind],
                    self.v[t, ind],
                    self.theta[t, ind],
                    self.dCons[t, ind],
                    self.FB1Cons[t, ind],
                    self.FB2Cons[t, ind]
                ) = self.stage_t_state_k_problem(t,k)
                if ind in self.hurricaneData.absorbing_states:
                    self.theta[t, ind].setAttr(GRB.Attr.UB, 0);
                    if self.inputParams.absorbing_option == 0:
                        for i in range(Ni):
                            for j in range(N0):
                                if j != i:
                                    self.f[t, ind][j,i].setAttr(GRB.Attr.UB, 0);

    # Train model: forward pass
    def FOSDDP_forward_pass_oneSP_iteration(self):
        k_init = self.inputParams.k_init;
        Ni = self.networkData.Ni;
        T = self.hurricaneData.T;
        absorbing_states = self.hurricaneData.absorbing_states;
        k_t = k_init-1;
        in_sample = [k_t];
        xval = np.zeros((Ni, T));
        thetaval = [0]*T;
        lb = 1e10;
        for t in range(T):
            if t > 0:
                k_t = MC_sample(in_sample[t - 1]+1, self.hurricaneData)
                in_sample.append(k_t)
                self.MSP_fa_update_RHS(k_t, t, xval)
            self.m[t, k_t].optimize()
            if self.m[t, k_t].status != GRB.OPTIMAL:
                print("Error in Forward Pass")
                print(f"Model in stage = {t} and state = {k_t}, in forward pass is {self.m[t, k_t].status}")
                sys.exit(0)
            else:
                for i in range(Ni):
                    xval[i, t] = self.x[t,k_t][i].x
                thetaval[t] = self.theta[t, k_t].x
                if t == 0:
                    lb = self.m[t, k_t].objVal
            if k_t in absorbing_states:
                break
        return xval, thetaval, lb, in_sample


    # Train model: backward pass
    def FOSDDP_backward_pass_oneSP_iteration(self, xval, thetaval, in_sample):
        T = self.hurricaneData.T;
        P_joint = self.hurricaneData.P_joint;
        Ni = self.networkData.Ni;
        K = self.hurricaneData.K;
        nodeLists = self.hurricaneData.nodeLists;
        absorbing_states = self.hurricaneData.absorbing_states;
        cutviolFlag = False;
        for t in range(len(in_sample)-1, 0, -1):
            # Solving all stage-t problems
            Q = [0] * K  #list for all the optimal values
            #list for all the dual multiplies of the first and second set of constraints
            pi1 = [[0] * Ni for _ in range(K)]
            pi2 = [[0] * Ni for _ in range(K)]
            sample_n = in_sample[t-1]; # the state observed at time t-1
            for k in range(len(nodeLists[t])):
                self.MSP_fa_update_RHS(nodeLists[t][k], t, xval)
                self.m[t, nodeLists[t][k]].optimize()
                if self.m[t, nodeLists[t][k]].status != GRB.OPTIMAL:
                    print("Error in Backward Pass")
                    print(f"Model in stage = {t} and state = {nodeLists[t][k]}, in backward pass is {self.m[t, nodeLists[t][k]].status}")
                    sys.exit(0)
                else:
                    Q[nodeLists[t][k]] = self.m[t, nodeLists[t][k]].objVal
                    for i in range(Ni):
                        pi1[nodeLists[t][k]][i] = self.FB1Cons[t,nodeLists[t][k]][i].pi
                        pi2[nodeLists[t][k]][i] = self.FB2Cons[t,nodeLists[t][k]][i].pi

            # Solving all stage-(t-1) problems and generate cuts/valid inequalities
            for n in range(len(nodeLists[t - 1])):
                if nodeLists[t - 1][n] not in absorbing_states:
                    Qvalue = 0
                    for k in range(len(nodeLists[t])):
                        if P_joint[nodeLists[t - 1][n]][nodeLists[t][k]] > self.hurricaneData.smallestTransProb:
                            Qvalue += Q[nodeLists[t][k]] * P_joint[nodeLists[t - 1][n]][nodeLists[t][k]]              					# check if cut is violated at the sample path encountered in the forward pass
                    
                    # check if cut is violated at the sample path encountered in the forward pass
                    if nodeLists[t - 1][n] == sample_n and (
                            (Qvalue - thetaval[t - 1]) / max(1e-10, abs(thetaval[t - 1])) > self.solveParams.cutviol
                            and abs(Qvalue - thetaval[t - 1]) > self.solveParams.cutviol
                    ):
                        cutviolFlag = True;
                    cutcoef = [0] * Ni
                    cutrhs_xval = 0
                    for k in range(len(nodeLists[t])):
                        if P_joint[nodeLists[t - 1][n]][nodeLists[t][k]] > self.hurricaneData.smallestTransProb:
                            for i in range(Ni):
                                tempval = (pi1[nodeLists[t][k]][i]+pi2[nodeLists[t][k]][i]) * P_joint[nodeLists[t - 1][n]][nodeLists[t][k]]
                                cutcoef[i] += tempval
                                cutrhs_xval += tempval * xval[i][t - 1]
                    self.m[t - 1, nodeLists[t - 1][n]].addConstr(
                        self.theta[t - 1, nodeLists[t - 1][n]] - gp.quicksum(
                            cutcoef[i] * self.x[t - 1, nodeLists[t - 1][n]][i] for i in range(Ni)
                        ) >= Qvalue - cutrhs_xval
                    )
        return cutviolFlag

    # Train model
    def train_models_offline(self):
        x_0 = self.networkData.x_0;
        T = self.hurricaneData.T;
        Ni = self.networkData.Ni;
        k_init = self.inputParams.k_init;
        # Set the RHS of the first_stage problem
        for i in range(Ni):
            self.FB1Cons[0,k_init-1][i].setAttr(GRB.Attr.RHS, x_0[i]);
            self.FB2Cons[0,k_init-1][i].setAttr(GRB.Attr.RHS, x_0[i]);

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
            xval, thetaval, lb, in_sample = self.FOSDDP_forward_pass_oneSP_iteration()
            LB.append(lb)
            # Termination check
            flag, Elapsed = termination_check(iter, relative_gap, LB, start, cutviol_iter, self.solveParams)
            if flag != 0:
                train_time = Elapsed
                break
            # Backward pass (if not terminated)
            cutviolFlag = self.FOSDDP_backward_pass_oneSP_iteration(xval, thetaval, in_sample)
            if cutviolFlag:
                cutviol_iter = 0
            else:
                cutviol_iter += 1
        return LB, train_time, iter

    # Update RHS of flow-balance and demand constraint
    def MSP_fa_update_RHS(self, k_t, t, xval):
        Ni = self.networkData.Ni;
        Nj = self.networkData.Nj;
        S = self.hurricaneData.states;
        T = self.hurricaneData.T;
        SCEN = self.networkData.SCEN;
        for i in range(Ni):
            self.FB1Cons[t,k_t][i].setAttr(GRB.Attr.RHS, xval[i,t - 1])
            self.FB2Cons[t,k_t][i].setAttr(GRB.Attr.RHS, xval[i,t - 1])
        for j in range(Nj):
            if S[k_t][2] == T and S[k_t][0] != 1:
                self.dCons[t,k_t][j].setAttr(GRB.Attr.RHS, SCEN[k_t][j]);
            else:
                self.dCons[t,k_t][j].setAttr(GRB.Attr.RHS, 0);

    # Evaluate model
    def FOSDDP_eval(self, osfname):
        Ni = self.networkData.Ni;
        Nj = self.networkData.Nj;
        N0 = self.networkData.N0;
        T = self.hurricaneData.T;
        nbOS = self.inputParams.nbOS;
        absorbing_states = self.hurricaneData.absorbing_states;

        self.define_models()

        LB, train_time, iter = self.train_models_offline()

        OS_paths = pd.read_csv(osfname).values
        objs_fa = np.zeros((nbOS, T))
        xval_fa = {}
        fval_fa = {}
        zval_fa = {}
        vval_fa = {}

        # key KPIs
        procurmnt_all = np.zeros((nbOS,T));
        procurmnt_amount = np.zeros(T); 
        procurmnt_percentage = np.zeros(T); 
        procurmnt_posExpect = np.zeros(T); 
        flow_amount = np.zeros(T);
        invAmount = np.zeros(nbOS);
        salvageAmount = np.zeros(nbOS);
        penaltyAmount = np.zeros(nbOS);

        start = time.time()
        for s in range(nbOS):
            xval = np.zeros((Ni, T))
            for t in range(T):
                k_t = OS_paths[s, t]-1
                if t > 1:
                    self.MSP_fa_update_RHS(k_t, t, xval)
                self.m[t,k_t].optimize()
                if self.m[t,k_t].status != GRB.OPTIMAL:
                    print(" in evaluation")
                    print(f"Model in stage = {t} and state = {k_t}, in forward pass is {self.m[t,k_t].status}")
                    exit(0)
                else:
                    objs_fa[s, t] = self.m[t,k_t].ObjVal - self.theta[t,k_t].x
                    xval_fa[s,t] = {}
                    for i in range(Ni):
                        xval[i, t] = self.x[t,k_t][i].x
                        xval_fa[s, t][i] = self.x[t,k_t][i].x
                    fval_fa[s, t] = {}
                    for i in range(N0):
                        for ii in range(Ni):
                            fval_fa[s, t][i,ii] = self.f[t,k_t][i,ii].x
                    zval_fa[s, t] = {}
                    for j in range(Nj):
                        zval_fa[s, t][j] = self.z[t,k_t][j].x
                    vval_fa[s, t] = {}
                    for i in range(Ni):
                        vval_fa[s, t][i] = self.v[t,k_t][i].x
                salvageAmount[s] += sum(vval_fa[s,t][i] for i in range(Ni));
                penaltyAmount[s] += sum(zval_fa[s,t][j] for j in range(Nj));
                procurmnt_amount[t] += sum(fval_fa[s,t][N0-1,i] for i in range(Ni));
                procurmnt_all[s,t] = sum(fval_fa[s,t][N0-1,i] for i in range(Ni));
                flow_amount[t] += sum(sum(fval_fa[s,t][i,ii] for i in range(Ni)) for ii in range(Ni));
                if k_t in absorbing_states:
                    invAmount[s] = sum(xval_fa[s,t-1][i] for i in range(Ni));
                    break
              
        fa_bar = np.mean(np.sum(objs_fa, axis=1))
        fa_std = np.std(np.sum(objs_fa, axis=1))
        fa_low = fa_bar - 1.96 * fa_std / np.sqrt(nbOS)
        fa_high = fa_bar + 1.96 * fa_std / np.sqrt(nbOS)
        CI = fa_bar-fa_low;
        print("FA...")
        print(f"μ ± 1.96*σ/√NS = {fa_bar} ± {CI}")
        test_time = time.time() - start

        # Now let's compute some KPIs
        for t in range(T):
            procurmnt_amount[t] = procurmnt_amount[t]/nbOS;
            flow_amount[t] = flow_amount[t]/nbOS;      
        for t in range(T):
            count = 0;
            totalPos = 0;
            for s in range(nbOS):
                if procurmnt_all[s,t] > 1e-2:
                    count += 1;
                    totalPos += procurmnt_all[s,t];
            procurmnt_percentage[t] = count*1.0/nbOS;
            if count > 0:
                procurmnt_posExpect[t] = totalPos*1.0/count;

        print("procurement amount = ", procurmnt_amount);
        print("flow amount = ", flow_amount);

        avgInvAmount = sum(invAmount[s] for s in range(nbOS))*1.0/nbOS;
        avgSalvageAmount = sum(salvageAmount[s] for s in range(nbOS))*1.0/nbOS;
        avgPenaltyAmount = sum(penaltyAmount[s] for s in range(nbOS))*1.0/nbOS;

        print("avgInvAmount = ", avgInvAmount);
        print("avgSalvageAmount = ", avgSalvageAmount);
        print("avgPenaltyAmount = ", avgPenaltyAmount);

        KPIvec = procurmnt_amount.tolist()+procurmnt_percentage.tolist()+procurmnt_posExpect.tolist()+flow_amount.tolist()+[avgInvAmount,avgSalvageAmount,avgPenaltyAmount]
        return [fa_bar, CI, train_time, test_time], KPIvec