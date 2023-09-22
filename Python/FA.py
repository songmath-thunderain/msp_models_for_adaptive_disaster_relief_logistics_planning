#MSP fully adaptive model functions 

#Define stage-t problem
def stage_t_state_k_problem(t):
    #######################
    # Define the model.
    m = Model('model')

    #######################
    # Define the variables.
    x = m.addVars(Ni, lb=0, ub=x_cap)
    f = m.addVars(N0, Ni, lb=0, ub=f_cap)
    y = m.addVars(Ni, Nj, lb=0)
    z = m.addVars(Nj, lb=0)
    v = m.addVars(Ni, lb=0)
    theta = m.addVar(lb = 0)

    #######################
    # Define the objective.
    m.setObjective(
        quicksum(cb[i, ii, t] * f[i, ii] for ii in range(1, Ni + 1) for i in range(1, N0 + 1)) +
        quicksum(ch[i, t] * x[i] for i in range(1, Ni + 1)) +
        quicksum(f[N0, i] for i in range(1, Ni + 1)) * h[t] +
        quicksum(ca[i, j, t] * y[i, j] for j in range(1, Nj + 1) for i in range(1, Ni + 1)) +
        quicksum(z[j] for j in range(1, Nj + 1)) * p +
        quicksum(v[i] for i in range(1, Ni + 1)) * q +
        theta,
        GRB.MINIMIZE
    )

    #######################
    # Define the constraints.
    FB1Cons = {}
    FB2Cons = {}
    dCons = {}

    for i in range(1, Ni + 1):
        if t == 1:
            FB1Cons[i] = m.addConstr(
                x[i] +
                quicksum(f[i, j] for j in range(1, Ni + 1) if j != i) -
                quicksum(f[j, i] for j in range(1, N0 + 1) if j != i) +
                quicksum(y[i, j] for j in range(1, Nj + 1)) +
                v[i] ==
                x_0[i]
            )
            FB2Cons[i] = m.addConstr(
                quicksum(f[i, j] for j in range(1, Ni + 1) if j != i) <= x_0[i]
            )
        else:
            FB1Cons[i] = m.addConstr(
                x[i] +
                quicksum(f[i, j] for j in range(1, Ni + 1) if j != i) -
                quicksum(f[j, i] for j in range(1, N0 + 1) if j != i) +
                quicksum(y[i, j] for j in range(1, Nj + 1)) +
                v[i] ==
                0
            )
            FB2Cons[i] = m.addConstr(
                quicksum(f[i, j] for j in range(1, Ni + 1) if j != i) <= 0
            )

    for j in range(1, Nj + 1):
        dCons[j] = m.addConstr(z[j] + quicksum(y[i, j] for i in range(1, Ni+1)) >= 0);
        
    return m, x, f, y, z, v, theta, dCons, FB1Cons, FB2Cons

###############################################################
###############################################################

def define_models():
    m, x, f, y, z, v, theta, dCons, FB1Cons, FB2Cons = {}, {}, {}, {}, {}, {}, {}, {}, {}, {}
    for t in range(1, T+1):
        if t == 1:
            m[t,k_init], x[t,k_init], f[t,k_init], y[t,k_init], z[t,k_init], v[t,k_init], theta[t,k_init], dCons[t,k_init], FB1Cons[t,k_init], FB2Cons[t,k_init] = stage_t_state_k_problem(t)
        else:
            for k in range(1, K+1):
                if k in absorbing_states:
                    continue
                else:
                    m[t,k], x[t,k], f[t,k], y[t,k], z[t,k], v[t,k], theta[t,k], dCons[t,k], FB1Cons[t,k], FB2Cons[t,k] = stage_t_state_k_problem(t)
    return m, x, f, y, z, v, theta, dCons, FB1Cons, FB2Cons


###############################################################
###############################################################
#Train model: forward pass
def FOSDDP_forward_pass_oneSP_iteration(lb, xval, thetaval):
    k_t = k_init.copy()
    in_sample = [k_t]  # what is the state in the first stage
    for t in range(1, T+1):
        # the state is known in the first stage; if not sample a new state k
        if t > 1:
            # sample a new state k
            k_t = MC_sample(in_sample[t-2])
            in_sample.append(k_t)
            # if k_t is absorbing no need to do any computation
            if k_t in absorbing_states:
                continue
            # update the RHS
            # MSP_fa_update_RHS(k_t,t,xval,rand(1:M)) # we do not have this second layer now [REVISION]
            MSP_fa_update_RHS(k_t, t, xval)
            
        # solve the model
        m_fa[t, k_t].optimize()
        
        # check the status
        status = m_fa[t, k_t].status
        if status != GRB.OPTIMAL:
            print("Error in Forward Pass")
            print("Model in stage =", t, "and state =", k_t, ", in forward pass is", status)
            exit(0)
        else:
            # collect values
            xval[:, t-1] = x_fa[t, k_t].x
            thetaval[t-1] = ϴ_fa[t, k_t].x
            if t == 1:
                lb = m_fa[t, k_t].ObjVal
    return xval, thetaval, lb, in_sample

###############################################################
###############################################################

#Train model: backward pass
def FOSDDP_backward_pass_oneSP_iteration(lb,xval,thetaval,in_sample):
    cutviolFlag = 0
    for t in range(T, 1, -1):
        # initialize
        Q = [0]*K # list for all the optimal values
        # list for all the dual multiplies of the first and second set of constraints
        pi1 = [[0]*Ni for _ in range(K)]
        pi2 = [[0]*Ni for _ in range(K)]
        sample_n = in_sample[t-2] # the states observed at time t-1
        for k in range(K):
            if k in absorbing_states:
                Q[k] = 0
                continue
            else:
                # Here we just update xval
                MSP_fa_update_RHS(k,t,xval)
                # solve the model
                m_fa[t,k].optimize()

                # check the status
                status = m_fa[t,k].status
                if status != GRB.OPTIMAL:
                    print("Error in Backward Pass")
                    print(f"Model in stage = {t} and state = {k}, in forward pass is {status}")
                    exit(0)
                else:
                    # collect values
                    Q[k] = m_fa[t,k].ObjVal
                    for i in range(Ni):
                        pi1[k][i] = FB1Cons_fa[t,k][i].pi
                        pi2[k][i] = FB2Cons_fa[t,k][i].pi

        for n in range(K):
            if n not in absorbing_states:
                if t-1 == 1 and n != k_init:
                    continue
                # what is the expected cost value
                Qvalue = sum(Q[k]*P_joint[n][k] for k in range(K))

                # check if cut is violated at the sample path encountered in the forward pass
                if n == sample_n and (Qvalue-thetaval[t-2])/max(1e-10,abs(thetaval[t-2])) > eps and abs(Qvalue-thetaval[t-2]) > eps:
                    cutviolFlag = 1

                # we are doing cut sharing so we will add the cut regardless
                m_fa[t-1,n].addConstr(ϴ_fa[t-2,n] - sum(sum((pi1[k][i]+pi2[k][i])*x_fa[t-2,n][i] for i in range(Ni))*P_joint[n][k] for k in range(K)) >= Qvalue - sum(sum((pi1[k][i]+pi2[k][i])*xval[i][t-2] for i in range(Ni))*P_joint[n][k] for k in range(K)))

    return cutviolFlag


###############################################################
###############################################################

#Train model
def train_models_offline():
    # set up model
    model = gp.Model()

    # set the RHS of the first_stage problem
    for i in range(Ni):
        FB1Cons_fa[1,k_init][i].setAttr(GRB.Attr.RHS, x_0[i])
        FB2Cons_fa[1,k_init][i].setAttr(GRB.Attr.RHS, x_0[i])

    # intialize stuff
    train_time = 0
    relative_gap = 1e10
    lb = 0
    LB = []
    xval = np.zeros((Ni,T))
    thetaval = np.zeros(T)
    iter = 0
    cutviol_iter = 0
    start = time()
    while True:
        iter += 1

        # forward pass
        xval, thetaval, lb, in_sample = FOSDDP_forward_pass_oneSP_iteration(lb, xval, thetaval)
        LB.append(lb)

        # termination check
        flag, Elapsed = termination_check(iter, relative_gap, LB, start, cutviol_iter)
        if flag != 0:
            train_time = Elapsed
            break

        # backward pass (if not terminated)
        cutviolFlag = FOSDDP_backward_pass_oneSP_iteration(lb, xval, thetaval, in_sample)
        if cutviolFlag == 1:
            cutviol_iter = 0
        else:
            cutviol_iter = cutviol_iter + 1

    return LB, train_time, iter

###############################################################
###############################################################

#evaluate model
# evaluate model
def FOSDDP_eval_offline():
    start = time.time()
    osfname = "JuliaData/OOS" + str(k_init) + ".csv"
    OS_paths = pd.read_csv(osfname, header=None).to_numpy() # read the out-of-sample file
    
    objs_fa = np.zeros((nbOS, T))
    
    xval_fa = [[None for _ in range(T)] for _ in range(nbOS)]
    fval_fa = [[None for _ in range(T)] for _ in range(nbOS)]
    yval_fa = [[None for _ in range(T)] for _ in range(nbOS)]
    zval_fa = [[None for _ in range(T)] for _ in range(nbOS)]
    vval_fa = [[None for _ in range(T)] for _ in range(nbOS)]
    
    procurmnt_amount = np.zeros(T)
    
    for s in range(nbOS):
        xval = np.zeros((Ni, T))
        for t in range(T):
            # the state is known in the first stage; if not sample a new state k
            k_t = int(OS_paths[s, t])
            # we do not have this second layer now [REVISION]
            # m = OS_M[s]; # realization from OS path corresponding to layer 2
            
            if k_t not in absorbing_states:
                if t > 1:
                    MSP_fa_update_RHS(k_t, t, xval)
                # solve the model
                m_fa_[t,k_t].optimize()
                # check the status
                status = m_fa_[t,k_t].status
                if status != GRB.OPTIMAL:
                    print(" in evaluation")
                    print("Model in stage =", t, " and state = ", k_t, ", in forward pass is ", status)
                    exit(0)
                else:
                    # collect values
                    xval_fa[s][t] = m_fa[t,k_t].getAttr('x', x_fa[t,k_t])
                    xval[:, t] = xval_fa[s][t]
                    fval_fa[s][t] = m_fa[t,k_t].getAttr('x', f_fa[t,k_t])
                    yval_fa[s][t] = m_fa[t,k_t].getAttr('x', y_fa[t,k_t])
                    zval_fa[s][t] = m_fa[t,k_t].getAttr('x', z_fa[t,k_t])
                    vval_fa[s][t] = m_fa[t,k_t].getAttr('x', v_fa[t,k_t])
                    objs_fa[s][t] = m_fa[t,k_t].objVal - m_fa[t,k_t].getAttr('x', theta_fa[t,k_t])
                    
                    procurmnt_amount[t] += np.sum(fval_fa[s][t][N0, :]) / nbOS
                            
        elapsed = time.time() - start
    fa_bar = np.mean(np.sum(objs_fa, axis=0))
    fa_std = np.std(np.sum(objs_fa, axis=0))
    fa_low = fa_bar - 1.96 * fa_std / np.sqrt(nbOS)
    fa_high = fa_bar + 1.96 * fa_std / np.sqrt(nbOS)
    print("FA...")
    print("mean = ", fa_bar, " interval = ", [fa_low, fa_high])
    elapsed = time() - start
    return objs_fa, fa_bar, fa_low, fa_high, elapsed

###############################################################
###############################################################

#update RHS of flow-balance and demand constraint
def MSP_fa_update_RHS(k_t, t, xval):
    for i in range(Ni):
        FB1Cons_fa[t,k_t][i].setAttr("rhs", xval[i,t-1])
        FB2Cons_fa[t,k_t][i].setAttr("rhs", xval[i,t-1])
    for j in range(Nj):
        if S[k_t][3] == Nc-1 and k_t not in absorbing_states:
            dCons_fa[t,k_t][j].setAttr("rhs", SCEN[k_t][j])
        else:
            dCons_fa[t,k_t][j].setAttr("rhs", 0)