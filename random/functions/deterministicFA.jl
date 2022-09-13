#MSP fully adaptive model functions, but assuming that the landfall time is deterministic

function non_terminal_stage_single_period_problem(t)
    #######################
    #Define the model.
    m = Model(optimizer_with_attributes(() -> Gurobi.Optimizer(GRB_ENV), "OutputFlag" => 0));

    #######################
    #Define the variables.
    @variables(m,
            begin
               0 <= x[i=1:Ni] <= x_cap[i];
               0 <= f[i=1:N0,ii=1:Ni] <= f_cap[i,ii]; 
               0 <= ϴ;
            end
          );

    #######################
    #Define the objective.
    @objective(m,
               Min,
               sum(sum(cb[i,ii,t]*f[i,ii] for ii=1:Ni) for i=1:N0)
              +sum(ch[i,t]*x[i] for i=1:Ni)
              +sum(f[N0,i] for i=1:Ni)*h[t]
              +ϴ
               );

    #######################
    #Define the constraints.
    
    #initialize two arrays to store the constraints which containt x_{t-1}
    FB1 = Array{Any,1}(undef,Ni);
    FB2 = Array{Any,1}(undef,Ni);
    
    #create the following constraints for every SP i
    #initialize the RHS to be zero for now, and will change it later 
    for i=1:Ni
        FB1[i] = @constraint(m, 
                             x[i]
                            +sum(f[i,j] for j=1:Ni if j != i)
                            -sum(f[j,i] for j=1:N0 if j != i)
                            == 0
                            );
        FB2[i] = @constraint(m, 
                            sum(f[i,j] for j=1:Ni if j != i)
                            <= 0
                            );
    end

    return m, x, f, ϴ, FB1, FB2
end

###############################################################
###############################################################

# This function defines the terminal stage for the MSP model with deterministic landfall
function terminal_stage_single_period_problem(t)
    #######################
    #Define the model.
    m = Model(optimizer_with_attributes(() -> Gurobi.Optimizer(GRB_ENV), "OutputFlag" => 0));

    #######################
    #Define the variables.
    @variables(m,
            begin
               0 <= x[i=1:Ni] <= x_cap[i];
               0 <= f[i=1:N0,ii=1:Ni] <= f_cap[i,ii];
               0 <= y[i=1:Ni,j=1:Nj];
               0 <= z[j=1:Nj];
            end
          );
    
    #######################
    #Define the objective.
    @objective(m,
               Min,
               sum(sum(cb[i,ii,t]*f[i,ii] for ii=1:Ni) for i=1:N0)
              +sum(ch[i,t]*x[i] for i=1:Ni)
              +sum(f[N0,i] for i=1:Ni)*h[t]
              +sum(sum(ca[i,j,t]*y[i,j] for j=1:Nj) for i=1:Ni)
              +sum(z[j] for j=1:Nj)*p
              +sum(x[i]-sum(y[i,j] for j=1:Nj) for i=1:Ni)*q
               );
    

    #######################
    #Define the constraints.
    
    #initialize two arrays to store the constraints which containt x_{t-1}
    FB1 = Array{Any,1}(undef,Ni);
    FB2 = Array{Any,1}(undef,Ni);
	dCons = Array{Any,1}(undef,Nj);
    
    #create the following constraints for every SP i
    #initialize the RHS to be zero for now, and will change it later 
    for i=1:Ni
        FB1[i] = @constraint(m, 
                             x[i]
                            +sum(f[i,j] for j=1:Ni if j != i)
                            -sum(f[j,i] for j=1:N0 if j != i)
                            == 0
                            );
        FB2[i] = @constraint(m, 
                            sum(f[i,j] for j=1:Ni if j != i)
                            <= 0
                            );
        @constraint(m,
                    sum(y[i,j] for j=1:Nj)
                    <=x[i]
                   );
    end
    for j=1:Nj
	   dCons[j] = @constraint(m, z[j]+sum(y[i,j] for i=1:Ni) == 0);
	end

#=
	for j=1:Nj
        # this constraint ensures that the flow sent to an SP is not more than the realized demand
        @constraint(m,
                    sum(y[i,j] for i=1:Ni)
                    <=d[j]
                   );
        # this constraint computes the unsatisfied demand
        @constraint(m,
                    d[j]-sum(y[i,j] for i=1:Ni)
                    <=z[j]
                   );
    end
=#

    # this constraint ensures that items cannot be shipped directly from the MDC to SP
    @constraint(m,
            sum(f[N0,i] for i=1:Ni) == 0
            );
    
    return m, x, f, y, z, FB1, FB2, dCons
end


###############################################################
###############################################################


# This function defines all the models for the MSP with deterministic landfall
function define_models(T)
    # first we initialize the list where we store all the models
    model = Array{Any,2}(undef,T,K); # list to store all the models for every stage and Markovian state
    x = Array{Any,2}(undef,T,K); # list to store all the x variables for every stage and Markovian state
    f = Array{Any,2}(undef,T,K); # ..............
    theta = Array{Any,2}(undef,T-1,K); # ..............
    y = Array{Any,1}(undef,K);
    z = Array{Any,1}(undef,K);
    FB1 = Array{Any,2}(undef,T,K);
    FB2 = Array{Any,2}(undef,T,K);
	dCons = Array{Any,1}(undef,K);
    
    #Then we define the a model for every stage and Markovian state state 
    for k=1:K, t=1:T
        if t < T
            # define use the function non_terminal_stage_single_period_problem
            model[t,k], x[t,k], f[t,k], theta[t,k], FB1[t,k], FB2[t,k] = non_terminal_stage_single_period_problem(t);
        else
            #use the function terminal_stage_single_period_problem
            model[t,k], x[t,k], f[t,k], y[k], z[k], FB1[t,k], FB2[t,k], dCons[k] = terminal_stage_single_period_problem(t);
        end
    end
    return model, x, f, theta, y, z, FB1, FB2, dCons
end

###############################################################
###############################################################

### TBD ###
#Train model: forward pass
function FOSDDP_forward_pass_oneSP_iteration(lb,xval,thetaval)
    k_t = copy(k_init);
    in_sample = [k_t]; #what is the state in the first stage
    for t=1:T
        #the state is known in the first stage; if not sample a new state k 
        if t>1
            #sample a new state k
            k_t = MC_sample(in_sample[t-1]);
            push!(in_sample,k_t);
            # if k_t is absorbing no need to do any computation
            if k_t in absorbing_states
                continue; 
            end
            #update the RHS
			#MSP_fa_update_RHS(k_t,t,xval,rand(1:M)); # we do not have this second layer now [REVISION]
			MSP_fa_update_RHS(k_t,t,xval);
        end
            
        #solve the model
        optimize!(m_fa[t,k_t]);

        #check the status 
        status = termination_status(m_fa[t,k_t]);
        if status != MOI.OPTIMAL
            println("Error in Forward Pass");
            println("Model in stage =", t, " and state = ", k_t, ", in forward pass is ", status);
            exit(0);
        else
            #collect values
            xval[:,t] = value.(x_fa[t,k_t]);
            thetaval[t] = value(ϴ_fa[t,k_t]);
            if t == 1
                lb = objective_value(m_fa[t,k_t]);
            end
        end
    end
    return xval, thetaval, lb, in_sample
end

###############################################################
###############################################################

#Train model: backward pass: new version, without two-layer uncertainty
function FOSDDP_backward_pass_oneSP_iteration(lb,xval,thetaval,in_sample)
    cutviolFlag = 0;
    for t=T:-1:2
        #initialize
        Q = zeros(K); #list for all the optimal values
         #list for all the dual multiplies of the first and second set of constraints
        pi1 = zeros(K,Ni); 
		pi2 = zeros(K,Ni); 
        sample_n = in_sample[t-1]; #the states observed at time t-1
        for k=1:K
            if k in absorbing_states
                Q[k] = 0;
                continue
            else
                # Here we just update xval
                MSP_fa_update_RHS(k,t,xval);
                #solve the model
                optimize!(m_fa[t,k]);

                #check the status 
                status = termination_status(m_fa[t,k]);
                if status != MOI.OPTIMAL
                    println("Error in Backward Pass");
                    println("Model in stage =", t, " and state = ", k, ", in forward pass is ", status);
                    exit(0);
                else
                    #collect values
                    Q[k] = objective_value(m_fa[t,k]);
                    for i=1:Ni
                        pi1[k,i] = dual(FB1Cons_fa[t,k][i]);
                        pi2[k,i] = dual(FB2Cons_fa[t,k][i]);
                    end
                end                 
            end
        end

        for n = 1:K
            if  n ∉ absorbing_states
                if t-1 == 1 && n != k_init
                    continue
                end
                #what is the expected cost value 
                Qvalue = sum(Q[k]*P_joint[n,k]  for k=1:K);

                # check if cut is violated at the sample path encountered in the forward pass
                if n == sample_n && (Qvalue-thetaval[t-1])/max(1e-10,abs(thetaval[t-1])) > ϵ
                    cutviolFlag = 1;
                end

                # we are doing cut sharing so we will add the cut regardless
                @constraint(m_fa[t-1,n],
                ϴ_fa[t-1,n]
                -sum(sum((pi1[k,i]+pi2[k,i])*x_fa[t-1,n][i] for i=1:Ni)*P_joint[n,k] for k=1:K)
                >=
                Qvalue-sum(sum((pi1[k,i]+pi2[k,i])*xval[i,t-1] for i=1:Ni)*P_joint[n,k] for k=1:K)
                );
            end
        end
    end
    return cutviolFlag;
end


###############################################################
###############################################################

#Train model
function train_models_offline()
    #set the RHS of the first_stage problem
    for i=1:Ni
        set_normalized_rhs(FB1Cons_fa[1,k_init][i], x_0[i]);
        set_normalized_rhs(FB2Cons_fa[1,k_init][i], x_0[i]);
    end

    #intialize stuff
    train_time = 0;
    relative_gap = 1e10;
    lb=0;
    LB = [];
    xval = zeros(Ni,T);
    thetaval = zeros(T);
    iter = 0;
    cutviol_iter = 0;
    start=time();
    while true
        iter+=1;
        #forward pass
        xval, thetaval, lb, in_sample = FOSDDP_forward_pass_oneSP_iteration(lb,xval,thetaval);
        push!(LB,lb);
		#println("LB = ", lb)
        #termination check
        flag, Elapsed = termination_check(iter,relative_gap,LB,start,cutviol_iter);
        if flag != 0
            train_time = Elapsed;
            break;
        end

        #backward pass (if not terminated)
        cutviolFlag = FOSDDP_backward_pass_oneSP_iteration(lb,xval,thetaval,in_sample);
        if cutviolFlag == 1
            cutviol_iter = 0;
        else
            cutviol_iter = cutviol_iter + 1;
        end
    end
    return LB, train_time, iter
end

###############################################################
###############################################################

#evaluate model
function FOSDDP_eval_offline()
    start=time();
    OS_paths = Matrix(CSV.read("./data/OOS.csv",DataFrame)); #read the out-of-sample file
	# we do not have this second layer now [REVISION]
	#OS_M = Matrix(CSV.read("./data/inOOS.csv",DataFrame))[:,1] #read the second layer OOS
    objs_fa = zeros(nbOS,T);
    
    xval_fa = Array{Any,2}(undef,nbOS,T); fval_fa = Array{Any,2}(undef,nbOS,T);
    yval_fa = Array{Any,2}(undef,nbOS,T); zval_fa = Array{Any,2}(undef,nbOS,T); vval_fa = Array{Any,2}(undef,nbOS,T);

    procurmnt_amount = zeros(T); 
    
    for s=1:nbOS
        xval = zeros(Ni,T);
        for t=1:T
            #the state is known in the first stage; if not sample a new state k 
            k_t = OS_paths[s,t];
            # we do not have this second layer now [REVISION]
			#m = OS_M[s]; # realization from OS path corresponding to layer 2

            if k_t ∉ absorbing_states
                if t > 1
                    MSP_fa_update_RHS(k_t,t,xval);
                end
                #solve the model
                optimize!(m_fa[t,k_t]);

                #check the status 
                status = termination_status(m_fa[t,k_t]);
                if status != MOI.OPTIMAL
                    println(" in evaluation");
                    println("Model in stage =", t, " and state = ", k_t, ", in forward pass is ", status);
                    exit(0);
                else
                    #collect values
                    xval_fa[s,t] = value.(x_fa[t,k_t]); xval[:,t] = xval_fa[s,t];
                    fval_fa[s,t] = value.(f_fa[t,k_t]);                 
                    yval_fa[s,t] = value.(y_fa[t,k_t]);
                    zval_fa[s,t] = value.(z_fa[t,k_t]);
                    vval_fa[s,t] = value.(v_fa[t,k_t]);
                    objs_fa[s,t] = objective_value(m_fa[t,k_t])- value(ϴ_fa[t,k_t]);
                    
                    procurmnt_amount[t] += (sum(fval_fa[s,t][N0,i] for i=1:Ni))/nbOS;
                end
            end
        end        
    end
    fa_bar = mean(sum(objs_fa[:,t] for t=1:T));
    fa_std = std(sum(objs_fa[:,t] for t=1:T));
    fa_low = fa_bar-1.96*fa_std/sqrt(nbOS);
    fa_high = fa_bar+1.96*fa_std/sqrt(nbOS);
	println("FA...");
    println("μ ± 1.96*σ/√NS = ", fa_bar, " ± ", [fa_low,fa_high]);
    elapsed = time() - start;
    vals = [xval_fa, fval_fa, yval_fa, zval_fa, vval_fa];
    return objs_fa, fa_bar, fa_low, fa_high, elapsed#, vals
end


###############################################################
###############################################################

#update RHS of flow-balance and demand constraint
# we do not have this second layer now [REVISION]
function MSP_fa_update_RHS(k_t,t,xval)
    for i=1:Ni        
        set_normalized_rhs(FB1Cons_fa[t,k_t][i], xval[i,t-1]);
        set_normalized_rhs(FB2Cons_fa[t,k_t][i], xval[i,t-1]);
    end 
    for j=1:Nj
        if S[k_t][3] == Nc-1 && k_t ∉ absorbing_states
            set_normalized_rhs(dCons_fa[t,k_t][j], SCEN[k_t][j]);
        else
            set_normalized_rhs(dCons_fa[t,k_t][j], 0);
        end
    end
end

