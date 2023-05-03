#MSP fully adaptive model functions 

#Define stage-t problem
function stage_t_state_k_problem(t)
 
    #######################
    #Define the model.
    m = Model(optimizer_with_attributes(() -> Gurobi.Optimizer(GRB_ENV), "OutputFlag" => 0));

    #######################
    #Define the variables.
    @variables(m,
            begin
               0 <= x[i=1:Ni] <= x_cap[i]
               0 <= f[i=1:N0,ii=1:Ni]
               0 <= y[i=1:Ni,j=1:Nj]
               0 <= z[j=1:Nj]
               0 <= v[i=1:Ni]
               0 <= ϴ
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
                    +sum(v[i] for i=1:Ni)*q
                    +ϴ
               );
    
    #######################
    #Define the constraints.
    FB1Cons = Dict(); #a dictonary to store flow-balance constraints 1
    FB2Cons = Dict(); #a dictonary to store  flow-balance constraints 2
    dCons = Dict(); #a dictonary to store all the demand constraint

    for i=1:Ni
        if t == 1
            FB1Cons[i] = @constraint(m, 
                                        x[i]
                                        +sum(f[i,j] for j=1:Ni if j != i)
                                        -sum(f[j,i] for j=1:N0 if j != i)
                                        +sum(y[i,j] for j=1:Nj)
                                        +v[i]             
                                        == x_0[i]
                                    );
            FB2Cons[i] = @constraint(m, sum(f[i,j] for j=1:Ni if j != i) <= x_0[i]);
        else
            FB1Cons[i] = @constraint(m, 
                                        x[i]
                                        +sum(f[i,j] for j=1:Ni if j != i)
                                        -sum(f[j,i] for j=1:N0 if j != i)
                                        +sum(y[i,j] for j=1:Nj) 
                                        +v[i]                
                                        == 0
                                    );
            FB2Cons[i] = @constraint(m, sum(f[i,j] for j=1:Ni if j != i) <= 0);  
        end
    end
    for j=1:Nj
       dCons[j] = @constraint(m, z[j]+sum(y[i,j] for i=1:Ni) >= 0);
    end
    
    return m, x, f, y, z, v, ϴ, dCons, FB1Cons, FB2Cons
end

###############################################################
###############################################################

#Define the model
function define_models()
    m = Dict(); x = Dict(); f = Dict(); y = Dict(); z = Dict(); v = Dict(); 
    ϴ = Dict(); dCons = Dict(); FB1Cons = Dict(); FB2Cons = Dict(); 
    for t=1:T
		#if t == 1
		#    m[t,k_init], x[t,k_init], f[t,k_init], y[t,k_init], z[t,k_init], v[t,k_init],
		#    ϴ[t,k_init], dCons[t,k_init], FB1Cons[t,k_init], FB2Cons[t,k_init] = stage_t_state_k_problem(t);
		#else
			for k=1:length(nodeLists[t])
				#if k in absorbing_states
				#	continue 
				#else
				ind = nodeLists[t][k];
                m[t,ind], x[t,ind], f[t,ind], y[t,ind], z[t,ind], v[t,ind], ϴ[t,ind], dCons[t,ind], FB1Cons[t,ind], FB2Cons[t,ind] = stage_t_state_k_problem(t);
				if ind in absorbing_states
					@constraint(m[t,ind], ϴ[t,ind] == 0);
					if absorbing_option == 0
						# not allowing MDC/SP operation at the absorbing state
						for i=1:Ni
							for j=1:Ni
								if j!= i
									@constraint(m[t,ind], f[t,ind][i,j] == 0);
								end
							end
							for j=1:N0
								if j!= i
									@constraint(m[t,ind], f[t,ind][j,i] == 0);
								end
							end
						end
					end
				end
            	#end
			end			
		#end
    end    
    return m, x, f, y, z, v, ϴ, dCons, FB1Cons, FB2Cons
end

###############################################################
###############################################################

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
			#if k_t in absorbing_states
			#    continue; 
			#end
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
		# if k_t is absorbing, terminate the forward pass
		if k_t in absorbing_states
			break;
		end
    end
    return xval, thetaval, lb, in_sample
end

###############################################################
###############################################################

#Train model: backward pass: new version, without two-layer uncertainty
function FOSDDP_backward_pass_oneSP_iteration(lb,xval,thetaval,in_sample)
	# This implements the full cut-sharing scheme
    cutviolFlag = 0;
    for t=length(in_sample):-1:2
		########################################
		# Solving all stage-t problems
		########################################
        #initialize
        Q = zeros(K); #list for all the optimal values
         #list for all the dual multiplies of the first and second set of constraints
        pi1 = zeros(K,Ni); 
		pi2 = zeros(K,Ni); 
        sample_n = in_sample[t-1]; #the states observed at time t-1
        for k=1:length(nodeLists[t])
			#if k in absorbing_states
			#    Q[k] = 0;
			#    continue
			#else
                # Here we just update xval
                MSP_fa_update_RHS(nodeLists[t][k],t,xval);
                #solve the model
                optimize!(m_fa[t,nodeLists[t][k]]);

                #check the status 
                status = termination_status(m_fa[t,nodeLists[t][k]]);
                if status != MOI.OPTIMAL
                    println("Error in Backward Pass");
                    println("Model in stage =", t, " and state = ", nodeLists[t][k], ", in forward pass is ", status);
                    exit(0);
                else
                    #collect values
                    Q[nodeLists[t][k]] = objective_value(m_fa[t,nodeLists[t][k]]);
                    for i=1:Ni
                        pi1[nodeLists[t][k],i] = shadow_price(FB1Cons_fa[t,nodeLists[t][k]][i]);
                        pi2[nodeLists[t][k],i] = shadow_price(FB2Cons_fa[t,nodeLists[t][k]][i]);
                    end
                end                 
			#end
        end
		########################################
		# Solving all stage-(t-1) problems and generate cuts/valid inequalities
		########################################
        for n = 1:length(nodeLists[t-1])
            if  nodeLists[t-1][n] ∉ absorbing_states
                #what is the expected cost value 
                Qvalue = 0;
				#tempProb = 0;
				for k=1:length(nodeLists[t])
					if P_joint[nodeLists[t-1][n],nodeLists[t][k]] > smallestTransProb
						Qvalue = Qvalue + Q[nodeLists[t][k]]*P_joint[nodeLists[t-1][n],nodeLists[t][k]];
						#tempProb = tempProb + P_joint[nodeLists[t-1][n],nodeLists[t][k]];
					end
				end
				#if abs(tempProb-1) > 1e-7
				#	@printf("tempProb = %f\n", tempProb);
				#	@printf("node = %d, \n", nodeLists[t-1][n]);
				#	println("S[node] = ", S[nodeLists[t-1][n]]);
				#	println("P_joint[", nodeLists[t-1][n], "] = ", P_joint[nodeLists[t-1][n],:]);
				#	exit(0);
				#end	
                # check if cut is violated at the sample path encountered in the forward pass
                if nodeLists[t-1][n] == sample_n && (Qvalue-thetaval[t-1])/max(1e-10,abs(thetaval[t-1])) > ϵ && abs(Qvalue-thetaval[t-1]) > ϵ
                    cutviolFlag = 1;
                end

                # we are doing cut sharing so we will add the cut regardless
				cutcoef = zeros(Ni);
				cutrhs_xval = 0;
				for k=1:length(nodeLists[t])
					if P_joint[nodeLists[t-1][n],nodeLists[t][k]] > smallestTransProb
						for i=1:Ni
							tempval = (pi1[nodeLists[t][k],i]+pi2[nodeLists[t][k],i])*P_joint[nodeLists[t-1][n],nodeLists[t][k]];
							cutcoef[i] = cutcoef[i] + tempval;
							cutrhs_xval = cutrhs_xval + tempval*xval[i,t-1];
						end
					end
				end
                @constraint(m_fa[t-1,nodeLists[t-1][n]],
                ϴ_fa[t-1,nodeLists[t-1][n]]-sum(cutcoef[i]*x_fa[t-1,nodeLists[t-1][n]][i] for i=1:Ni) >= Qvalue-cutrhs_xval);
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
	osfname = "./data/OOS"*string(k_init)*".csv";
    OS_paths = Matrix(CSV.read(osfname,DataFrame)); #read the out-of-sample file
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
			if k_t in absorbing_states
				break;
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
        if S[k_t][3] == Nc && S[k_t][1] != 1
            set_normalized_rhs(dCons_fa[t,k_t][j], SCEN[k_t][j]);
        else
            set_normalized_rhs(dCons_fa[t,k_t][j], 0);
        end
    end
end

