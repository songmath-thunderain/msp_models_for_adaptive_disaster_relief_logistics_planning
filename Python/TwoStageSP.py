import time
import gurobipy as gp
from gurobipy import GRB
import pandas as pd
import numpy as np
from misc import *
import sys
import copy

class TwoStageSP:
    def __init__(self,inputParams,solveParams,hurricaneData,networkData):
        self.inputParams = inputParams;
        self.solveParams = solveParams;
        self.hurricaneData = hurricaneData;
        self.networkData = networkData;

    # Define the terminal-stage problem: only used when absorbing_option = 1, i.e., MDC/SP operation is allowed to occur
    def terminal_model(self, t_roll, x_init):
        # Define the model
        m = gp.Model()
        # Data instantiation
        Ni = self.networkData.Ni;
        N0 = self.networkData.N0;
        Nj = self.networkData.Nj;
        ca = self.networkData.ca;
        cb = self.networkData.cb;
        ch = self.networkData.ch;
        cp = self.networkData.cp;
        p = self.networkData.p;
        q = self.networkData.q;
        x_cap = self.networkData.x_cap;

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
    def RH_2SSP_first_stage(self, t_roll, k_t, x_init):
        # Note that the static 2SSP corresponds to the case when t_roll = 0
        # Data instantiation
        Ni = self.networkData.Ni;
        N0 = self.networkData.N0;
        Nj = self.networkData.Nj;
        ca = self.networkData.ca;
        cb = self.networkData.cb;
        ch = self.networkData.ch;
        cp = self.networkData.cp;
        p = self.networkData.p;
        q = self.networkData.q;
        x_cap = self.networkData.x_cap;
        T = self.hurricaneData.T;
        absorbing_option = self.inputParams.absorbing_option;
        nodeScenList = self.hurricaneData.nodeScenList;
        nodeScenWeights = self.hurricaneData.nodeScenWeights;

        nbstages1 = T - t_roll
        if absorbing_option == 0:
            # If no MDC/SP operation is allowed in the absorbing state, do not plan for stage T
            nbstages1 = T - t_roll - 1
        
        # Define the model
        m = gp.Model()
        nbScens = len(nodeScenList[(t_roll,k_t)])
        pbScens = nodeScenWeights[(t_roll,k_t)]
        
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
            theta[n] = m.addVar(lb=-1e8)
        
        # Define the objective
        m.setObjective(
                gp.quicksum(
                    gp.quicksum(cb[i,ii,t+t_roll] * f[i, ii, t] for ii in range(Ni) for i in range(N0))
                    for t in range(nbstages1))
            + gp.quicksum(
                gp.quicksum(ch[i,t_roll+t] * x[i, t] for i in range(Ni))
                for t in range(nbstages1))
            + gp.quicksum(
                gp.quicksum(cp[t_roll+t] * f[N0-1, i, t] for i in range(Ni))
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
        m.update();
        m.setParam("OutputFlag", 0);
        return m, x, f, theta

    # Define second-stage scenario subproblem
    def RH_2SSP_second_stage(self):
        # Data instantiation
        Ni = self.networkData.Ni;
        N0 = self.networkData.N0;
        Nj = self.networkData.Nj;
        ca = self.networkData.ca;
        cb = self.networkData.cb;
        ch = self.networkData.ch;
        cp = self.networkData.cp;
        p = self.networkData.p;
        q = self.networkData.q;
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
        m.update();
        m.setParam("OutputFlag", 0);
        return m, y, xCons, dCons, rCons


    # Define the two-stage SP models: master problem and subproblem
    def RH_2SSP_define_models(self, t_roll, k_t, x_init):
        # Define first stage (master problem) model
        self.master, self.x, self.f, self.theta = self.RH_2SSP_first_stage(t_roll, k_t, x_init)
        
        # Define second stage (subproblem) optimality model
        self.subproblem, self.y2, self.xCons, self.dCons, self.rCons = self.RH_2SSP_second_stage()
        

    # Solves the two-stage SP model
    def RH_2SSP_solve_roll(self, k_t, t_roll):
        Ni = self.networkData.Ni;
        N0 = self.networkData.N0;
        T = self.hurricaneData.T;
        absorbing_option = self.inputParams.absorbing_option;
        nodeScenList = self.hurricaneData.nodeScenList;
        nodeScenWeights = self.hurricaneData.nodeScenWeights;

        nbstages1 = T - t_roll
        if absorbing_option == 0:
            # If no MDC/SP operation is allowed in the absorbing state, do not plan for stage T since we know for sure that all states are absorbing
            nbstages1 = T - t_roll - 1
        nbScens = len(nodeScenList[(t_roll,k_t)])
        pbScens = nodeScenWeights[(t_roll,k_t)]
        LB = -1e10
        UB = 1e10
        thetaval = [0] * nbScens
        xval = [[0] * nbstages1 for _ in range(Ni)]
        fval = [[[0] * nbstages1 for _ in range(Ni)] for _ in range(N0)]
        solveIter = 0
        flag = 1
        
        while (UB - LB) * 1.0 / max(1e-10, abs(LB)) > self.solveParams.cutviol and abs(UB - LB) > self.solveParams.cutviol and solveIter < 100:  # WARNING: Temporary fix here!
            # Solve first stage
            flag = 0
            solveIter += 1
            LB, xval, fval, thetaval = self.solve_first_stage(nbstages1, nbScens)
            firstCost = LB - sum(thetaval[n] * pbScens[n] for n in range(nbScens))
            
            # Solve second stage
            flag, Qbar = self.solve_second_stage(k_t, t_roll, xval, fval, thetaval)
            if flag != -1:
                UB = min(firstCost + Qbar, UB)
        
        if solveIter == 100:
            print("# iterations is maxed out!")
            print("LB =", LB)
            print("UB =", UB)
        
        return LB, UB, xval, fval, thetaval

    # Solves the first-stage problem
    def solve_first_stage(self, nbstages1, nbScens):
        Ni = self.networkData.Ni;
        N0 = self.networkData.N0;
        self.master.optimize()
        status_master = self.master.status
        if status_master != GRB.OPTIMAL:
            print("Master problem status is:", status_master, "=== Oops! :/")
            sys.exit(0)
        else:
            LB = self.master.objVal
            xval = [[self.x[i, t].x for t in range(nbstages1)] for i in range(Ni)] # referred to as xval[i][t]
            fval = [[[self.f[i, j, t].x for t in range(nbstages1)] for j in range(Ni)] for i in range(N0)]
            thetaval = [self.theta[n].x for n in range(nbScens)]
            return LB, xval, fval, thetaval

    # Solves the second-stage problem
    def solve_second_stage(self, k_t, t_roll, xval, fval, thetaval):
        Ni = self.networkData.Ni;
        N0 = self.networkData.N0;
        T = self.hurricaneData.T;
        cb = self.networkData.cb;
        ch = self.networkData.ch;
        cp = self.networkData.cp;
        absorbing_option = self.inputParams.absorbing_option;
        nodeScenList = self.hurricaneData.nodeScenList;
        nodeScenWeights = self.hurricaneData.nodeScenWeights;
        
        nbScens = len(nodeScenList[(t_roll,k_t)])
        pbScens = nodeScenWeights[(t_roll,k_t)]
        flag = 0
        nbstages1 = T - t_roll
        
        if absorbing_option == 0:
            nbstages1 = T - t_roll - 1
        
        Q = [0] * nbScens
        pi1 = [[] for _ in range(nbScens)]
        pi2 = [[] for _ in range(nbScens)]
        pi3 = [0] * nbScens
        Qbar = 0
        
        for n in range(nbScens):
            absorbingT = nodeScenList[(t_roll,k_t)][n][0]
            absorbingState = nodeScenList[(t_roll,k_t)][n][1]
            self.RH_2SSP_update_RHS(absorbingT, absorbingState, xval, fval, t_roll)
            Q[n], pi1[n], pi2[n], pi3[n], flag = self.solve_scen_subproblem()
            
            if flag == -1:
                print("Subproblem status is infeasible?!")
                exit(0)
        
        Qbar = sum(Q[n] * pbScens[n] for n in range(nbScens))
        
        # Cut generation: multi-cut version
        for n in range(nbScens):
            #print("Q[%d] = %f, (%d, %d)" %(n, Q[n], nodeScenList[(t_roll,k_t)][n][0], nodeScenList[(t_roll,k_t)][n][1]))
            if (Q[n] - thetaval[n]) / max(1e-10, abs(Q[n])) > self.solveParams.cutviol and Q[n] - thetaval[n] > self.solveParams.cutviol:
                #print("pi1[n] = ", pi1[n]);
                #print("pi3[n] = %f" % pi3[n])
                tt = nodeScenList[(t_roll,k_t)][n][0]
                # tt is the terminal stage
                if absorbing_option == 0:
                    reimbursement = -sum(
                            sum(sum(
                                cb[i,ii,t_roll + t] * fval[i][ii][t]
                                for ii in range(Ni)
                            ) for i in range(N0))
                        + sum(ch[i,t_roll + t] * xval[i][t] for i in range(Ni))
                        + sum(fval[N0-1][i][t] for i in range(Ni)) * cp[t_roll + t]
                        for t in range(tt - t_roll, nbstages1));
                    #print("reimbursement[%d] = %f" %(n,reimbursement));
                    self.master.addConstr(
                        self.theta[n] - sum(pi1[n][i] * self.x[i, tt - t_roll - 1] for i in range(Ni))
                        + pi3[n] * (-sum(
                            sum(sum(
                                cb[i,ii,t_roll + t] * self.f[i,ii,t]
                                for ii in range(Ni)
                            ) for i in range(N0))
                        + sum(ch[i,t_roll + t] * self.x[i, t] for i in range(Ni))
                        + sum(self.f[N0-1,i,t] for i in range(Ni)) * cp[t_roll + t]
                        for t in range(tt - t_roll, nbstages1)
                        )) >= Q[n] - sum(pi1[n][i] * xval[i][tt - t_roll -1] for i in range(Ni)) + pi3[n]*reimbursement
                    )
                else:
                    reimbursement = -sum(
                            sum(sum(
                                cb[i,ii,t_roll + t] * fval[i][ii][t]
                                for ii in range(Ni)
                            ) for i in range(N0))
                        + sum(ch[i,t_roll + t] * xval[i][t] for i in range(Ni))
                        + sum(fval[N0-1][i][t] for i in range(Ni)) * cp[t_roll + t]
                        for t in range(tt + 1 - t_roll, nbstages1));
                    #print("reimbursement[%d] = %f" %(n,reimbursement));
                    self.master.addConstr(
                        self.theta[n] - sum(pi1[n][i] * self.x[i, tt - t_roll] for i in range(Ni))
                        + pi3[n] * (-sum(
                            sum(sum(
                                cb[i,ii,t_roll + t] * self.f[i,ii,t]
                                for ii in range(Ni)
                            ) for i in range(N0))
                        + sum(ch[i,t_roll + t] * self.x[i, t] for i in range(Ni))
                        + sum(self.f[N0-1,i,t] for i in range(Ni)) * cp[t_roll + t]
                        for t in range(tt + 1 - t_roll, nbstages1)
                        )) >= Q[n] - sum(pi1[n][i] * xval[i][tt - t_roll] for i in range(Ni)) + pi3[n]*reimbursement
                    )
                
                flag = 1
        
        return flag, Qbar

    # Solves scenario subproblem of the second stage
    def solve_scen_subproblem(self):
        flag = 0;
        Ni = self.networkData.Ni;
        Nj = self.networkData.Nj;
        self.subproblem.optimize()
        status_subproblem = self.subproblem.status
        
        Qtemp = 0
        pi1temp = [0] * Ni
        pi2temp = [0] * Nj
        pi3temp = 0
        
        if status_subproblem != GRB.OPTIMAL:
            if status_subproblem == GRB.INFEASIBLE:
                flag = -1
        else:
            Qtemp = self.subproblem.objVal
            pi1temp = [self.xCons[i].pi for i in range(Ni)]
            pi2temp = [self.dCons[j].pi for j in range(Nj)]
            pi3temp = self.rCons.pi
        
        return Qtemp, pi1temp, pi2temp, pi3temp, flag

    # Updates the RHS of the 2nd-stage constraints and objective coefficients
    def RH_2SSP_update_RHS(self, absorbingT, k_t, xval, fval, t_roll):
        Ni = self.networkData.Ni;
        Nj = self.networkData.Nj;
        N0 = self.networkData.N0;
        SCEN = self.networkData.SCEN;
        T = self.hurricaneData.T;
        cb = self.networkData.cb;
        ca = self.networkData.ca;
        ch = self.networkData.ch;
        cp = self.networkData.cp;
        absorbing_option = self.inputParams.absorbing_option;
        nbstages1 = T - t_roll
        
        if absorbing_option == 0:
            # If no MDC/SP operation is allowed in the absorbing state, do not plan for stage T since we know for sure that all states are absorbing
            nbstages1 = T - t_roll - 1
        
        for i in range(Ni):
            if absorbing_option == 0:
                self.xCons[i].setAttr(GRB.Attr.RHS, xval[i][absorbingT - t_roll - 1])
            else:
                self.xCons[i].setAttr(GRB.Attr.RHS, xval[i][absorbingT - t_roll])

        for j in range(Nj):
            if self.hurricaneData.states[k_t][0] != 1:
                self.dCons[j].setAttr(GRB.Attr.RHS, SCEN[k_t][j]);
            else:
                self.dCons[j].setAttr(GRB.Attr.RHS, 0);
        
        if absorbingT == (T-1):
            # Plan exactly until the landfall time -- no reimbursement occurred!
            self.rCons.setAttr(GRB.Attr.RHS, 0);
        else:
            updatedRHS = 0
            
            if absorbing_option == 0:
                # reimburse the operational cost starting from the terminal stage, since the terminal stage does not allow operation
                updatedRHS = -sum((
                    sum(sum(cb[i,ii,t_roll + t] * fval[i][ii][t] for ii in range(Ni)) for i in range(N0))
                    + sum(ch[i,t_roll + t] * xval[i][t] for i in range(Ni)) + sum(fval[N0-1][i][t] for i in range(Ni)) * cp[t_roll + t])
                    for t in range(absorbingT-t_roll,nbstages1))
            else:
                # reimburse the operational cost if they occur after the terminal stage: starting from stage (τ+1)-t_roll
                updatedRHS = -sum((
                    sum(sum(cb[i,ii,t_roll + t] * fval[i][ii][t] for ii in range(Ni)) for i in range(N0))
                    + sum(ch[i,t_roll + t] * xval[i][t] for i in range(Ni)) + sum(fval[N0-1][i][t] for i in range(Ni)) * cp[t_roll + t])
                    for t in range(absorbingT-t_roll+1,nbstages1))
            self.rCons.setAttr(GRB.Attr.RHS, updatedRHS);
        # Also need to update the coefficients of y[i,j] variables in the 2nd stage
        for i in range(Ni):
            for j in range(Nj):
                self.y2[i,j].setAttr(GRB.Attr.Obj, ca[i,j,absorbingT]);


    def static_2SSP_eval(self,osfname):
        T = self.hurricaneData.T;
        absorbing_option = self.inputParams.absorbing_option;
        nbOS = self.inputParams.nbOS;
        dissipate_option = self.inputParams.dissipate_option;
        nodeScenList = self.hurricaneData.nodeScenList;
        nodeScenWeights = self.hurricaneData.nodeScenWeights;
        S = self.hurricaneData.states;
        x_0 = self.networkData.x_0;

        s = 0
        t_roll = 0
        x_init = x_0

        # Define the model
        OS_paths = pd.read_csv(osfname).values 
        self.RH_2SSP_define_models(t_roll, OS_paths[s, t_roll]-1, x_init)

        # Solve the model
        start_time = time.time()
        LB, UB, xval, fval, thetaval = self.RH_2SSP_solve_roll(OS_paths[s, t_roll]-1, t_roll)
        timeTrain = time.time() - start_time

        nbScens = len(nodeScenList[t_roll, OS_paths[s, t_roll]-1])

        f1cost = LB - sum(thetaval[n] * nodeScenWeights[t_roll, OS_paths[s, t_roll]-1][n] for n in range(0, nbScens))

        print("training LB =", LB)
        print("training UB =", UB)

        # Evaluate the model
        start_time = time.time()

        objs = [f1cost] * nbOS
        Q = np.zeros(nbOS)
        pi1 = [None] * nbOS
        pi2 = [None] * nbOS
        pi3 = np.zeros(nbOS)

        for s in range(nbOS):
            absorbingT = -1
            if dissipate_option == 1:
                absorbingT = list(OS_paths[s, 0:T]).index(next((x for x in OS_paths[s, 0:T] if S[x-1][2] == T or S[x-1][0] == 1), None))
            else:
                absorbingT = list(OS_paths[s, 0:T]).index(next((x for x in OS_paths[s, 0:T] if S[x-1][2] == T), None))

            self.RH_2SSP_update_RHS(absorbingT, OS_paths[s, absorbingT]-1, xval, fval, t_roll)
            
            Q[s], pi1[s], pi2[s], pi3[s], flag = self.solve_scen_subproblem()
            objs[s] = objs[s] + Q[s]

        st2SSP_bar = np.mean(objs)
        st2SSP_std = np.std(objs)
        st2SSP_low = st2SSP_bar - 1.96 * st2SSP_std / np.sqrt(nbOS)
        st2SSP_high = st2SSP_bar + 1.96 * st2SSP_std / np.sqrt(nbOS)
        CI = 1.96 * st2SSP_std / np.sqrt(nbOS)

        print("static 2SSP....")
        print("μ ± 1.96*σ/√NS =", st2SSP_bar, "±", CI)
        timeTest = time.time() - start_time
        return [objs, st2SSP_bar, st2SSP_low, st2SSP_high, timeTrain, timeTest]

    def RH_2SSP_eval(self, osfname):
        T = self.hurricaneData.T;
        absorbing_option = self.inputParams.absorbing_option;
        nbOS = self.inputParams.nbOS;
        dissipate_option = self.inputParams.dissipate_option;
        nodeScenList = self.hurricaneData.nodeScenList;
        nodeScenWeights = self.hurricaneData.nodeScenWeights;
        S = self.hurricaneData.states;
        x_0 = self.networkData.x_0;
        SCEN = self.networkData.SCEN;
        Ni = self.networkData.Ni;
        Nj = self.networkData.Nj;
        N0 = self.networkData.N0;
        cb = self.networkData.cb;
        ca = self.networkData.ca;
        ch = self.networkData.ch;
        cp = self.networkData.cp;
        p = self.networkData.p;
        q = self.networkData.q;

        s = 0
        t_roll = 0
        x_init = x_0

        start_time = time.time()

        OS_paths = pd.read_csv(osfname).values 
        # Define and solve the model
        self.RH_2SSP_define_models(t_roll, OS_paths[s, t_roll]-1, x_init)
        LB_1, UB_1, xval_1, fval_1, thetaval_1 = self.RH_2SSP_solve_roll(OS_paths[s, t_roll]-1, t_roll)

        objs_RH2SSP = np.zeros((nbOS, T))

        for s in range(nbOS):
            for i in range(N0):
                for ii in range(Ni):
                    objs_RH2SSP[s, 0] += cb[i, ii, 0] * fval_1[i][ii][0]

            for i in range(Ni):
                objs_RH2SSP[s, 0] += (ch[i, 0] * xval_1[i][0] + cp[0] * fval_1[N0 - 1][i][0])

        for s in range(nbOS):
            absorbingT = -1
            if dissipate_option == 1:
                absorbingT = list(OS_paths[s, 0:T]).index(next((x for x in OS_paths[s, 0:T] if S[x-1][2] == T or S[x-1][0] == 1), None))
            else:
                absorbingT = list(OS_paths[s, 0:T]).index(next((x for x in OS_paths[s, 0:T] if S[x-1][2] == T), None))

            #x_init = [item[0] for item in xval_1]
            x_init = np.array(xval_1)[:,0]
            if dissipate_option == 1 and S[OS_paths[s, absorbingT]-1][0] == 1:
                # This rolling procedure will go all the way until the hurricane gets into the absorbing state of dissipating 
                for t_roll in range(1, absorbingT):
                    # roll up to t = absorbingT-1
                    self.RH_2SSP_define_models(t_roll, OS_paths[s, t_roll]-1, x_init)
                    LB_Roll, UB_Roll, xval_Roll, fval_Roll, thetaval_Roll = self.RH_2SSP_solve_roll(OS_paths[s, t_roll]-1, t_roll)
                    
                    #implement xₜ, pay cost, and pass xₜ to new t+1.
                    #x_init = [item[0] for item in xval_Roll]
                    x_init = np.array(xval_Roll)[:,0]
                    objs_RH2SSP[s, t_roll] = 0

                    for i in range(N0):
                        for ii in range(Ni):
                            objs_RH2SSP[s, t_roll] += cb[i, ii, t_roll] * fval_Roll[i][ii][0]

                    for i in range(Ni):
                        objs_RH2SSP[s, t_roll] += (ch[i, t_roll] * xval_Roll[i][0] + cp[t_roll] * fval_Roll[N0 - 1][i][0])

                    if t_roll == (absorbingT - 1):
                        # Now we get the realization, do the recourse now and finish the rolling procedure
                        t_roll += 1
                        objs_RH2SSP[s, t_roll] = 0
                        # Just salvage all the x_init

                        # Regardless of what the absorbing_option is, we won't do anything but to salvage everything
                        for i in range(Ni):
                            objs_RH2SSP[s, t_roll] += xval_Roll[i][0] * q
            else:
                # This rolling procedure will stop at t = absorbingT
                for t_roll in range(1, absorbingT):
                    # roll up to t = absorbingT-1
                    self.RH_2SSP_define_models(t_roll, OS_paths[s, t_roll]-1, x_init)
                    LB_Roll, UB_Roll, xval_Roll, fval_Roll, thetaval_Roll = self.RH_2SSP_solve_roll(OS_paths[s, t_roll]-1, t_roll)
                    
                    #x_init = [item[0] for item in xval_Roll]
                    x_init = np.array(xval_Roll)[:,0]
                    objs_RH2SSP[s, t_roll] = 0

                    for i in range(N0):
                        for ii in range(Ni):
                            objs_RH2SSP[s, t_roll] += cb[i, ii, t_roll] * fval_Roll[i][ii][0]

                    for i in range(Ni):
                        objs_RH2SSP[s, t_roll] += (ch[i, t_roll] * xval_Roll[i][0] + cp[t_roll] * fval_Roll[N0 - 1][i][0])

                    if t_roll == (absorbingT - 1):
                        t_roll += 1
                        objs_RH2SSP[s, t_roll] = 0

                        # Approach #2: based on xvals_Roll[:,0], optimize all the operations together with full information
                        if absorbing_option == 0:
                            for i in range(Ni):
                                self.xCons[i].setAttr(GRB.Attr.RHS, xval_Roll[i][0])
                            k_t = OS_paths[s, t_roll]-1
                            for j in range(Nj):
                                if S[k_t][0] != 1:
                                    self.dCons[j].setAttr(GRB.Attr.RHS, SCEN[k_t][j])
                                else:
                                    self.dCons[j].setAttr(GRB.Attr.RHS, 0)

                            self.rCons.setAttr(GRB.Attr.RHS, 0) #This is 0 since the cost in the future has not been paid yet -- this is rolling horizon
                            # Also need to update the coefficients of y2[i,j] variables in the 2nd stage
                            for i in range(Ni):
                                for j in range(Nj):
                                    self.y2[i, j].setAttr(GRB.Attr.Obj, ca[i, j, absorbingT])

                            self.subproblem.optimize()

                            if self.subproblem.status != GRB.OPTIMAL:
                                print("status =", self.subproblem.status)
                                exit(0)
                            else:
                                objs_RH2SSP[s, t_roll] += self.subproblem.objVal
                        else:
                            # Solve a terminal stage problem, just as the FA/MSP version
                            m_term, x_term, f_term, y_term, z_term, v_term, dCons_term = self.terminal_model(t_roll, x_init)
                            k_t = OS_paths[s, t_roll]-1

                            for j in range(Nj):
                                if S[k_t][0] != 1:
                                    dCons_term[j].setAttr(GRB.Attr.RHS, SCEN[k_t][j])
                                else:
                                    dCons_term[j].setAttr(GRB.Attr.RHS, 0)

                            m_term.optimize()

                            if m_term.status != GRB.OPTIMAL:
                                print("status =", m_term.status)
                                sys.exit(0)
                            else:
                                objs_RH2SSP[s, t_roll] += m_term.objVal

        elapsed_RH2SSP = time.time() - start_time

        RH2SSP_bar = np.mean(np.sum(objs_RH2SSP[:, :T], axis=1))
        RH2SSP_std = np.std(np.sum(objs_RH2SSP[:, :T], axis=1))
        RH2SSP_low = RH2SSP_bar - 1.96 * RH2SSP_std / np.sqrt(nbOS)
        RH2SSP_high = RH2SSP_bar + 1.96 * RH2SSP_std / np.sqrt(nbOS)
        CI = 1.96 * RH2SSP_std / np.sqrt(nbOS)
        print("RH 2SSP....")
        print("μ ± 1.96*σ/√NS =", RH2SSP_bar, "±", CI)
        return [objs_RH2SSP, RH2SSP_bar, RH2SSP_low, RH2SSP_high, elapsed_RH2SSP]

    def WS_eval(self, osfname):
        T = self.hurricaneData.T;
        absorbing_option = self.inputParams.absorbing_option;
        nbOS = self.inputParams.nbOS;
        absorbing_states = self.hurricaneData.absorbing_states;
        dissipate_option = self.inputParams.dissipate_option;
        nodeLists = self.hurricaneData.nodeLists;
        nodeScenList = self.hurricaneData.nodeScenList;
        nodeScenWeights = self.hurricaneData.nodeScenWeights;
        smallestTransProb = self.hurricaneData.smallestTransProb; 
        P_joint = self.hurricaneData.P_joint; 
        S = self.hurricaneData.states;
        x_0 = self.networkData.x_0;
        SCEN = self.networkData.SCEN;
        Ni = self.networkData.Ni;
        Nj = self.networkData.Nj;
        N0 = self.networkData.N0;
        cb = self.networkData.cb;
        ca = self.networkData.ca;
        ch = self.networkData.ch;
        cp = self.networkData.cp;
        p = self.networkData.p;
        q = self.networkData.q;
        
        start_time = time.time()

        # Initialization
        decisionNodes = copy.deepcopy(nodeLists);
        objvalNodes = [[] for _ in range(T)]
        solutionNodes = {}

        for t in range(T):
            objvalNodes[t] = [0.0] * len(nodeLists[t])

        for t_roll in range(T):
            for k in range(len(nodeLists[t_roll])):
                if nodeLists[t_roll][k] not in absorbing_states:
                    # transient state: solving a 2SSP
                    x_init = x_0
                    decisionNodes[t_roll][k] = 1
                    self.RH_2SSP_define_models(t_roll, nodeLists[t_roll][k], x_init)
                    LB_Roll, UB_Roll, xval_Roll, fval_Roll, thetaval_Roll = self.RH_2SSP_solve_roll(nodeLists[t_roll][k], t_roll)
                    objvalNodes[t_roll][k] = UB_Roll # objGo: expected objval if implementing a two-stage plan now. This is also only temporary, need to do a round of backward cleanup
                    solutionNodes[(t_roll, k)] = [xval_Roll, fval_Roll]
                else:
                    # absorbing state, we can determine pretty easily if it is Go vs. No-Go
                    if S[nodeLists[t_roll][k]][0] == 1 and dissipate_option == 1:
                        objvalNodes[t_roll][k] = 0
                        decisionNodes[t_roll][k] = 0
                    else:
                        # Hurricane makes landfall with intensity, need to decide to do nothing or do some last-minute operation (although perhaps costly)
                        costNoGo = p * np.sum(SCEN[nodeLists[t_roll][k]])
                        if absorbing_option == 0:
                            objvalNodes[t_roll][k] = costNoGo
                            decisionNodes[t_roll][k] = 0
                        else:
                            #define terminal stage optimality model
                            m_term, x_term, f_term, y_term, z_term, v_term, dCons_term = self.terminal_model(t_roll, x_0)
                            for j in range(Nj):
                                dCons_term[j].setAttr(GRB.Attr.RHS, SCEN[nodeLists[t_roll][k]][j])

                            m_term.optimize()

                            if m_term.status != GRB.OPTIMAL:
                                print("status_subproblem =", m_term.status)
                                exit(0)
                            else:
                                costGo = m_term.ObjVal

                            if costNoGo < costGo:
                                objvalNodes[t_roll][k] = costNoGo
                                decisionNodes[t_roll][k] = 0
                            else:
                                objvalNodes[t_roll][k] = costGo
                                decisionNodes[t_roll][k] = 1

        # Now let's do a round of backward correction!
        for t in range(T - 2, -1, -1):
            for k in range(len(nodeLists[t])):
                if nodeLists[t][k] not in absorbing_states:
                    costNoGo = 0
                    for kk in range(len(nodeLists[t + 1])):
                        if P_joint[nodeLists[t][k]][nodeLists[t + 1][kk]] > smallestTransProb:
                            costNoGo += P_joint[nodeLists[t][k]][nodeLists[t + 1][kk]] * objvalNodes[t + 1][kk]
                    if costNoGo < objvalNodes[t][k]:
                        decisionNodes[t][k] = 0
                        objvalNodes[t][k] = costNoGo

        train_time = time.time() - start_time
        start_time = time.time()

        # Start evaluating policies on the decision tree
        print("Construction is done....Now we do evaluation...")

        OS_paths = pd.read_csv(osfname).values 
        objs_OOS = np.zeros(nbOS)

        for s in range(nbOS):
            for t in range(T):
                ind = list(nodeLists[t]).index(next((x for x in nodeLists[t] if x == (OS_paths[s,t]-1)), None))
                if (OS_paths[s, t]-1) in absorbing_states:
                    # if absorbing, just take whatever that is the best, which has been computed above
                    objs_OOS[s] = objvalNodes[t][ind]
                    #print("absorbed! obj = ", objs_OOS[s], "\n");
                    break
                else:
                    if decisionNodes[t][ind] == 0:
                        # Decision is No-go!
                        continue
                    else:
                        # Decision is Go!
                        absorbingT = -1
                        if dissipate_option == 1:
                            absorbingT = list(OS_paths[s, 0:T]).index(next((x for x in OS_paths[s, 0:T] if S[x-1][2] == T or S[x-1][0] == 1), None))
                        else:
                            absorbingT = list(OS_paths[s, 0:T]).index(next((x for x in OS_paths[s, 0:T] if S[x-1][2] == T), None))

                        #define second stage (subproblem) optimality model
                        self.subproblem, self.y2, self.xCons, self.dCons, self.rCons = self.RH_2SSP_second_stage()

                        if absorbing_option == 0:
                            for i in range(Ni):
                                self.xCons[i].setAttr(GRB.Attr.RHS, solutionNodes[(t, ind)][0][i][absorbingT - t - 1])
                        else:
                            for i in range(Ni):
                                self.xCons[i].setAttr(GRB.Attr.RHS, solutionNodes[(t, ind)][0][i][absorbingT - t])

                        for j in range(Nj):
                            self.dCons[j].setAttr(GRB.Attr.RHS, SCEN[OS_paths[s,absorbingT]-1][j]);
                        
                        for i in range(Ni):
                            for j in range(Nj):
                                self.y2[i, j].setAttr(GRB.Attr.Obj, ca[i, j, absorbingT])

                        self.subproblem.optimize()

                        if subproblem.status != GRB.OPTIMAL:
                            print("status_subproblem =", subproblem.status)
                            exit(0)
                        else:
                            objs_OOS[s] = subproblem.ObjVal
                            #print("first obj = ", objs_OOS[s]);
                            if absorbing_option == 0:
                                for tt in range(absorbingT - t):
                                    objs_OOS[s] += sum(sum(
                                        cb[i, ii, t + tt] * solutionNodes[(t, ind)][1][i][ii][tt] for ii in range(Ni)) for i in range(N0)
                                    ) + sum(ch[i, t + tt] * solutionNodes[(t, ind)][0][i][tt] for i in range(Ni)) + sum(
                                        solutionNodes[(t, ind)][1][N0-1][i][tt] for i in range(Ni)
                                    ) * cp[t + tt]
                            else:
                                for tt in range(absorbingT + 1 - t):
                                    objs_OOS[s] += sum(sum(
                                        cb[i, ii, t + tt] * solutionNodes[(t, ind)][1][i][ii][tt] for ii in range(Ni)) for i in range(N0)
                                    ) + sum(ch[i, t + tt] * solutionNodes[(t, ind)][0][i][tt] for i in range(Ni)) + sum(
                                        solutionNodes[(t, ind)][1][N0-1][i][tt] for i in range(Ni)
                                    ) * cp[t + tt]
                        #print("Go! obj = ", objs_OOS[s], "\n");
                        break

        test_time = time.time() - start_time

        WS_bar = np.mean(objs_OOS)
        WS_std = np.std(objs_OOS)
        WS_low = WS_bar - 1.96 * WS_std / np.sqrt(nbOS)
        WS_high = WS_bar + 1.96 * WS_std / np.sqrt(nbOS)
        print("WS....")
        print("μ ± 1.96*σ/√NS =", WS_bar, "±", [WS_low, WS_high])