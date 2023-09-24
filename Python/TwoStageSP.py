#Define first-stage master problem
def RH_2SSP_first_stage(t_roll, x_init):
    # Note that the static 2SSP corresponds to the case when t_roll = 1
    nbstages1 = T - t_roll + 1
    m = Model('2SSP')
    
    # Define the variables, note that t is a relative index (to the current roll) going from 1 to nbstages1; 
    # from relative index to absolute index: t-> t_roll-1+t
    x = m.addVars(Ni, nbstages1, lb=0, ub=x_cap)
    f = m.addVars(N0, Ni, nbstages1, lb=0, ub=f_cap)
    theta = m.addVars(nbscen, lb=-1e8)
    
    # Define the objective.
    m.setObjective(
        quicksum(quicksum(quicksum(cb[i, ii, t_roll-1+t]*f[i, ii, t] for ii in range(Ni)) for i in range(N0)) for t in range(nbstages1))
        + quicksum(quicksum(ch[i, t_roll-1+t]*x[i, t] for i in range(Ni)) for t in range(nbstages1))
        + quicksum(quicksum(h[t_roll-1+t]*f[N0, i, t] for i in range(Ni)) for t in range(nbstages1))
        + quicksum(theta[n]*1.0/nbscen for n in range(nbscen)), GRB.MINIMIZE)
    
    # Define the constraints.
    for t in range(nbstages1):
        for i in range(Ni):
            if t == 0:
                m.addConstr(x[i, t] + quicksum(f[i, j, t] for j in range(Ni) if j != i)
                            - quicksum(f[j, i, t] for j in range(N0) if j != i) == x_init[i])
                m.addConstr(quicksum(f[i, j, t] for j in range(Ni) if j != i) <= x_init[i])
            else:
                m.addConstr(x[i, t-1] - x[i, t] - quicksum(f[i, j, t] for j in range(Ni) if j != i)
                            + quicksum(f[j, i, t] for j in range(N0) if j != i) == 0)
                m.addConstr(quicksum(f[i, j, t] for j in range(Ni) if j != i) <= x[i, t-1])
    
    m.Params.OutputFlag = 0
    m.update()
    
    return m, x, f, theta

###############################################################
###############################################################

#Define second-stage scenario supbproblem
def RH_2SSP_second_stage():
    #######################
    m = Model()
    #######################
    #Define the variables.
    y = m.addVars(Ni, Nj, lb=0)
    z = m.addVars(Nj, lb=0)
    v = m.addVars(Ni, lb=0)
    reimbursement = m.addVar(lb=-gp.GRB.INFINITY, ub=0)

    #######################
    #Define the objective.
    m.setObjective(
        quicksum(ca[i,j,1]*y[i,j] for j in range(Nj) for i in range(Ni))
        + quicksum(z[j] for j in range(Nj))*p
        + quicksum(v[i] for i in range(Ni))*q
        + reimbursement,
        GRB.MINIMIZE
    )

    #######################
    #Define the constraints.
    xCons = {} #a dictonary to store all the inventory constraints
    dCons = {} #a dictonary to store all the demand constraints

    for i in range(Ni):
        xCons[i] = m.addConstr(quicksum(y[i,j] for j in range(Nj)) + v[i] == 0)

    for j in range(Nj):
        dCons[j] = m.addConstr(quicksum(y[i,j] for i in range(Ni)) + z[j] == 0)

    rCons = m.addConstr(reimbursement == 0)

    return m, y, xCons, dCons, rCons

###############################################################
###############################################################

#defines the two-stage SP models: master problem and subproblem
def RH_2SSP_define_models(t_roll,x_init):
    #define first stage (master problem) model  
    master, x, f, theta = RH_2SSP_first_stage(t_roll,x_init);
    
    #define second stage (subproblem) optimality model
    subproblem, y2, xCons, dCons, rCons = RH_2SSP_second_stage();
    
    return master, x, f, theta, subproblem, y2, xCons, dCons, rCons
end


###############################################################
###############################################################

#initialize parameter for two-stage model
def initialize(s, t_roll):
    # s is the index for the sample path in the out-of-sample test 
    nbstages1 = T - t_roll + 1
    LB = -1e10
    UB = 1e10
    thetaval = np.zeros(nbscen, dtype=np.float64)
    xval = np.zeros((Ni, nbstages1), dtype=np.float64)
    fval = np.zeros((N0, Ni, nbstages1), dtype=np.float64)
    # Do sampling to create (in-sample) scenarios 
    osfname = "JuliaData/OOS" + str(k_init) + ".csv"
    allscen = pd.read_csv(osfname).to_numpy()  # read the out-of-sample file
    # In the first roll, always choose the nbscen scenarios out of the total of 10000 OOS once every 10000/nbscen 
    scen = allscen[range(0, 10000, 10000 // nbscen), :T]  # note that this is just an initialization for scen

    # when t_roll = 1, i.e., the first roll, it will always be the evenly selected scenarios from allscen

    # In later rolls, create in-sample scenarios by sampling
    if t_roll > 1:
        for n in range(nbscen):
            for t in range(1, t_roll):
                scen[n, t] = allscen[s, t]

            for t in range(t_roll, T):
                scen[n, t] = MC_sample(scen[n, t - 1])

    qprob = np.full(nbscen, 1 / nbscen, dtype=np.float64)

    return LB, UB, xval, fval, thetaval, scen, qprob

###############################################################
###############################################################

#solves the two-stage SP model
def RH_2SSP_solve_roll(s, t_roll, master, subproblem, x, f, theta, y, xCons, dCons, rCons):
    nbstages1 = T-t_roll+1
    LB, UB, xval, fval, θval, scen, qprob = initialize(s,t_roll)
    solveIter = 0
    while (UB-LB)*1.0/max(1e-10,abs(LB)) > eps and abs(UB-LB) > eps and solveIter < 100: # WARNING: Temporary fix here!
        # solve first stage
        solveIter += 1
        LB, xval, fval, thetaval = solve_first_stage(LB,xval,fval,thetaval,master,x,f,theta)
        firstCost = LB - sum(thetaval[n]*qprob[n] for n in range(1,nbscen+1))
        if t_roll < T:
            # solve second stage
            flag, Qbar = solve_second_stage(t_roll,xval,fval,thetaval,scen,qprob,master,subproblem,x,f,theta,y,xCons,dCons,rCons)
            if flag != -1:
                UB = min(firstCost+Qbar,UB)
        else:
            print("exiting here")
            break
    if solveIter == 100:
        print("# iterations is maxed out!")
        print("LB = ", LB)
        print(", UB = ", UB)
return LB, UB, xval, fval, thetaval

###############################################################
###############################################################

#solves the first-stage problem
def solve_first_stage(LB, xval, fval, θval, master, x, f, theta):
    master.optimize() # solve the model
    status_master = master.status # check the status
    if status_master != GRB.OPTIMAL:
        print("Master problem status is: ", status_master, " === Oops! :/")
        return
    else:
        # update the values
        LB = master.objVal
        xval = x.X
        fval = f.X
        thetaval = theta.X
    return LB, xval, fval, thetaval

###############################################################
###############################################################

#solves the second-stage problem
def solve_second_stage(t_roll,xval,fval,thetaval,scen,qprob,master,subproblem,x,f,theta,y,xCons,dCons,rCons):
    flag = 0
    nbstages1 = T-t_roll+1
    Q = [0]*nbscen
    pi1 = [[]]*nbscen
    pi2 = [[]]*nbscen
    pi3 = [0]*nbscen
    Qbar = 0
    absorbingT = -1
    tau = None
    
    for n in range(nbscen):
        #identify the period when the hurricane makes landfall 
        tau = next((x for x in range(len(scen[n,:])) if S[x][3] == Nc-1 and x not in absorbing_states), None)

        if tau is None:
            absorbingT = next((x for x in range(len(scen[n,:])) if S[x][1] == 1), None)
            print("n = ", n)
            print(", tau = ", tau)
            print(", absorbingT = ", absorbingT)
            RH_2SSP_update_RHS(absorbingT,scen[n,absorbingT],subproblem,xCons,dCons,rCons,xval,fval,y,t_roll)
        else:
            absorbingT = -1
            RH_2SSP_update_RHS(tau,scen[n,tau],subproblem,xCons,dCons,rCons,xval,fval,y,t_roll)

        #solve the subproblem and store the dual information
        Q[n], pi1[n], pi2[n], pi3[n], flag = solve_scen_subproblem(subproblem,xCons,dCons,rCons)

        if flag == -1:
            print("subproblem status is infeasible?!")
            exit(0)

    Qbar = sum(Q[n]*qprob[n] for n in range(nbscen))

    # cut generation: multi-cut version
    for n in range(nbscen):
        if (Q[n]-thetaval[n])/max(1e-10,abs(Q[n])) > eps and abs(Q[n]-thetaval[n]) > eps:
            tau = next((x for x in range(len(scen[n,:])) if S[x][3] == Nc-1 and x not in absorbing_states), None)
            if tau is None:
                tt = next((x for x in range(len(scen[n,:])) if S[x][1] == 1), None)
            else:
                tt = tau
            master.addConstr(theta[n]-quicksum(pi1[n][i]*x[i,tt-t_roll+1] for i in range(Ni))-pi3[n]*(-quicksum((sum(sum(cb[i,ii,t_roll+t-1]*f[i,ii,t] for ii in range(Ni)) for i in range(N0))+sum(ch[i,t_roll+t-1]*x[i,t] for i in range(Ni))+sum(f[N0,i,t] for i in range(Ni))*h[t_roll+t-1]) for t in range(tt+2-t_roll,nbstages1)) >= Q[n]-quicksum(pi1[n][i]*xval[i,tt-t_roll+1] for i in range(Ni))-pi3[n]*(-quicksum((sum(sum(cb[i,ii,t_roll+t-1]*fval[i,ii,t] for ii in range(Ni)) for i in range(N0))+sum(ch[i,t_roll+t-1]*xval[i,t] for i in range(Ni))+sum(fval[N0,i,t] for i in range(Ni))*h[t_roll+t-1]) for t in range(tt+2-t_roll,nbstages1))
            flag = 1;

    return flag, Qbar

###############################################################
###############################################################
#solves scenario subproblem of the second stage
def solve_scen_subproblem(subproblem, xCons, dCons, rCons):
    flag = 0
    subproblem.optimize() #solve the model
    status_subproblem = subproblem.status #check the status 
    Qtemp = 0
    pi1temp = [0]*Ni
    pi2temp = [0]*Nj
    pi3temp = 0
    if status_subproblem != GRB.OPTIMAL:
        if status_subproblem == GRB.INFEASIBLE:
            flag = -1
    else:
        #update the values
        Qtemp = subproblem.objVal
        pi1temp = [xCons[i].pi for i in range(Ni)]
        pi2temp = [dCons[j].pi for j in range(Nj)]
        pi3temp = rCons.pi
    return Qtemp, pi1temp, pi2temp, pi3temp, flag

###############################################################
###############################################################

#updates the RHS of the 2nd-stage constraints and objective coefficients
def RH_2SSP_update_RHS(tau, k_t, subproblem, xCons, dCons, rCons, xval, fval, y, t_roll):
    # Note that the assumption coming into this function is that tau > t_roll if tau is not None
    nbstages1 = T - t_roll + 1
    for i in range(Ni):
        if tau is None:
            print("We shouldn't come here any more!")
            exit(0)
            # set_normalized_rhs(xCons[i],xval[i,T-t_roll+1]);
        else:
            set_normalized_rhs(xCons[i], xval[i, tau - t_roll + 1])
    for j in range(Nj):
        if k_t not in absorbing_states:
            set_normalized_rhs(dCons[j], SCEN[k_t][j])
        else:
            set_normalized_rhs(dCons[j], 0)

    if tau is None:
        # nothing to reimburse here
        print("We shouldn't come here any more!")
        exit(0)
        # set_normalized_rhs(rCons, 0)
    else:
        updatedRHS = -sum((sum(sum(cb[i, ii, t_roll + t - 1] * fval[i, ii, t] for ii in range(1, Ni + 1)) for i in range(1, N0 + 1))
                       + sum(ch[i, t_roll + t - 1] * xval[i, t] for i in range(1, Ni + 1))
                       + sum(fval[N0, i, t] for i in range(1, Ni + 1)) * h[t_roll + t - 1])
                      for t in range(tau + 2 - t_roll, nbstages1 + 1))
        set_normalized_rhs(rCons, updatedRHS)
        # Also need to update the coefficients of y[i,j] variables in the 2nd stage
        for i in range(1, Ni + 1):
            for j in range(1, Nj + 1):
                set_objective_coefficient(subproblem, y[i, j], ca[i, j, tau])
