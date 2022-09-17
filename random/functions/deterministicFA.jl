#MSP fully adaptive model functions, but assuming that the landfall time is deterministic

function non_terminal_stage_single_period_problem_FAD(t)
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

	if t == Tmin	
		# this constraint ensures that items cannot be shipped directly from the MDC to SP
  		@constraint(m, sum(f[N0,i] for i=1:Ni) == 0);
	end

    return m, x, f, ϴ, FB1, FB2
end

###############################################################
###############################################################

# This function defines the terminal stage T for the MSP model with deterministic landfall
function terminal_stage_single_period_problem_FAD()
	#######################
    #Define the model.
    m = Model(optimizer_with_attributes(() -> Gurobi.Optimizer(GRB_ENV), "OutputFlag" => 0));

    #######################
    #Define the variables.
    @variables(m,
            begin
               0 <= y[i=1:Ni,j=1:Nj];
               0 <= z[j=1:Nj];
			   0 <= v[i=1:Ni];
            end
          );
    
    #######################
    #Define the objective.
    @objective(m,
               Min,
              sum(sum(ca[i,j,Tmin]*y[i,j] for j=1:Nj) for i=1:Ni)
              +sum(z[j] for j=1:Nj)*p
              +sum(v[i] for i=1:Ni)*q
               );
    

    #######################
    #Define the constraints.
    
    #initialize two arrays to store the constraints which contains x_{T}
    FB = Array{Any,1}(undef,Ni);
	dCons = Array{Any,1}(undef,Nj);
    
    #create the following constraints for every SP i
    #initialize the RHS to be zero for now, and will change it later 
    for i=1:Ni
        FB[i] = @constraint(m, sum(y[i,j] for j=1:Nj)+v[i] == 0);
    end
    for j=1:Nj
	   	dCons[j] = @constraint(m, z[j]+sum(y[i,j] for i=1:Ni) == 0);
	end

    return m, FB, dCons, y
end


###############################################################
###############################################################

# This function defines all the models for the MSP with deterministic landfall
function define_models_FAD()
    # first we initialize the list where we store all the models
    model = Array{Any,2}(undef,T,K); # list to store all the models for every stage and Markovian state
    x = Array{Any,2}(undef,T,K); # list to store all the x variables for every stage and Markovian state
    f = Array{Any,2}(undef,T,K); # ..............
    theta = Array{Any,2}(undef,T,K); # ..............
    FB1 = Array{Any,2}(undef,T,K);
    FB2 = Array{Any,2}(undef,T,K);
    
    #Then we define the a model for every stage and Markovian state state 
    for k=1:K, t=1:Tmin
        # define use the function non_terminal_stage_single_period_problem
        model[t,k], x[t,k], f[t,k], theta[t,k], FB1[t,k], FB2[t,k] = non_terminal_stage_single_period_problem_FAD(t);
    end

	model_final = Array{Any,1}(undef,K);
	y_final = Array{Any,1}(undef,K);
    FB_final = Array{Any,1}(undef,K);
	dCons_final = Array{Any,1}(undef,K);

	#To accomodate deterministic landfall model for a random landfall time, need to define a "final" stage problem
	for k=1:K
		model_final[k], FB_final[k], dCons_final[k], y_final[k] = terminal_stage_single_period_problem_FAD();
	end

    return model, x, f, theta, FB1, FB2, model_final, FB_final, dCons_final, y_final
end

###############################################################
###############################################################

#Train model: forward pass, up to Tmin, we assume that there is no demand realization yet. If landfall does occur on Tmin, then there will be a deterministic realization, otherwise, there will be a random realization
function FOSDDP_forward_pass_oneSP_iteration_FAD(lb,xval,thetaval)
    k_t = copy(k_init);
    in_sample = [k_t]; #what is the state in the first stage
    for t=1:Tmin
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
			for i=1:Ni
				set_normalized_rhs(FB1Cons_fa[t,k_t][i], xval[i,t-1]);
      	  		set_normalized_rhs(FB2Cons_fa[t,k_t][i], xval[i,t-1]);
			end
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

#Train model: backward pass
function FOSDDP_backward_pass_oneSP_iteration_FAD(lb,xval,thetaval,in_sample)
    cutviolFlag = 0;
    for t=(Tmin+1):-1:2
        sample_n = in_sample[t-1]; #the states observed at time t-1
		if t <= Tmin
			# There is nothing related to the demand realization/landfall
			#initialize
			Q = zeros(K); #list for all the optimal values
			 #list for all the dual multiplies of the first and second set of constraints
			pi1 = zeros(K,Ni); 
			pi2 = zeros(K,Ni);
			for k=1:K
				if k in absorbing_states
					Q[k] = 0;
					continue
				else
					# Here we just update xval
					for i=1:Ni
						set_normalized_rhs(FB1Cons_fa[t,k][i], xval[i,t-1]);
      	 		 		set_normalized_rhs(FB2Cons_fa[t,k][i], xval[i,t-1]);
					end
					#solve the model
					optimize!(m_fa[t,k]);

					#check the status 
					status = termination_status(m_fa[t,k]);
					if status != MOI.OPTIMAL
						println("Error in Backward Pass");
						println("Model in stage =", t, " and state = ", k, ", in backward pass is ", status);
						exit(0);
					else
						#collect values
						Q[k] = objective_value(m_fa[t,k]);
						for i=1:Ni
							pi1[k,i] = shadow_price(FB1Cons_fa[t,k][i]);
							pi2[k,i] = shadow_price(FB2Cons_fa[t,k][i]);
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
		else
			# This is the tricky part, need to consider the demand realization/landfall
			if sample_n in absorbing_states
				continue
			else
				if S[sample_n][3] == Nc-1
					# made landfall -> deterministic realization
					for i = 1:Ni
						set_normalized_rhs(FB_final[sample_n][i], xval[i,t-1]);
			   		end
					for j = 1:Nj
						set_normalized_rhs(dCons_final[sample_n][j], SCEN[sample_n][j]);
					end
					for i=1:Ni
						for j=1:Nj
							set_objective_coefficient(model_final[sample_n], y_final[sample_n][i,j], ca[i,j,Tmin]);
						end
					end
					#solve the model
					optimize!(model_final[sample_n]);

					#check the status 
					status = termination_status(model_final[sample_n]);
					if status != MOI.OPTIMAL
						println("Error in Backward Pass, final step");
						println("Model in stage =", t, " and state = ", sample_n, ", in backward pass is ", status);
						exit(0);
					else
						#collect values
						lastQ = objective_value(model_final[sample_n]);
						lastpi = zeros(Ni); 
						for i=1:Ni
							lastpi[i] = shadow_price(FB_final[sample_n][i]);
						end
						if (lastQ-thetaval[t-1])/max(1e-10,abs(thetaval[t-1])) > ϵ
							@constraint(m_fa[t-1,sample_n],
								ϴ_fa[t-1,sample_n]
								-sum(lastpi[i]*x_fa[t-1,sample_n][i] for i=1:Ni)
								>=
								lastQ-sum(lastpi[i]*xval[i,t-1] for i=1:Ni)
								);
						end
					end
				else
					# has not made landfall yet -> random realization, go ahead and do random sampling just as we did for RH	
					scen = zeros(Int,nbscen,T);
					for n=1:nbscen
						for tt=1:Tmin
							scen[n,tt] = in_sample[tt];
						end
						for tt=(Tmin+1):T
							scen[n,tt] = MC_sample(scen[n,tt-1]);
						end
					end
					lastQ = zeros(nbscen); #list for all the optimal values
					lastpi = Array{Any,1}(undef,nbscen); #list for all the dual multiplies of the first set of constraints
					qprob = fill(1/nbscen,nbscen);
					absorbingT = -1;
					τ = nothing
					for n=1:nbscen
						#identify the period when the hurricane makes landfall 
						τ = findfirst(x -> S[x][3] == Nc-1 && x ∉ absorbing_states, scen[n,:]);
						
						#update the RHS
						if τ === nothing     
							absorbingT = findfirst(x -> S[x][1] == 1, scen[n,:]);
							for i = 1:Ni
								set_normalized_rhs(FB_final[sample_n][i], xval[i,t-1]);
							end
							for j = 1:Nj
								set_normalized_rhs(dCons_final[sample_n][j], 0);
							end
							for i=1:Ni
								for j=1:Nj
									set_objective_coefficient(model_final[sample_n], y_final[sample_n][i,j], ca[i,j,absorbingT]);
								end
							end
							for tt=(Tmin+1):absorbingT
								lastQ[n] += sum(ch[i,tt]*xval[i,t-1] for i=1:Ni);
							end
						else
							for i = 1:Ni
								set_normalized_rhs(FB_final[sample_n][i], xval[i,t-1]);
							end
							for j = 1:Nj
								set_normalized_rhs(dCons_final[sample_n][j], SCEN[scen[n,τ]][j]);
							end
							for i=1:Ni
								for j=1:Nj
									set_objective_coefficient(model_final[sample_n], y_final[sample_n][i,j], ca[i,j,τ]);
								end
							end
							for tt=(Tmin+1):τ 
								lastQ[n] += sum(ch[i,tt]*xval[i,t-1] for i=1:Ni);
							end
						end
						#solve the subproblem and store the dual information
						optimize!(model_final[sample_n]) #solve the model
						status_subproblem = termination_status(model_final[sample_n]); #check the status 
						pitemp = zeros(Ni);
						if status_subproblem != MOI.OPTIMAL
							println("Error in Backward Pass, final step");
							println("Model in stage =", t, " and state = ", sample_n, ", in backward pass is ", status);
							exit(0);
						else
							#update the values
							lastQ[n] += objective_value(model_final[sample_n]);
							#need to include the inventory cost
							for i=1:Ni
								pitemp[i] = shadow_price(FB_final[sample_n][i]);
							end
							lastpi[n] = pitemp;
						end
					end
					lastQbar = sum(lastQ[n]*qprob[n] for n=1:nbscen);
					
					if (lastQbar-thetaval[t-1])/max(1e-10,abs(thetaval[t-1])) > ϵ
						@constraint(m_fa[t-1,sample_n],
								ϴ_fa[t-1,sample_n]
								-sum(qprob[n]*sum(lastpi[n][i]*x_fa[t-1,sample_n][i] for i=1:Ni) for n=1:nbscen)
								>=
								lastQbar-sum(qprob[n]*sum(lastpi[n][i]*xval[i,t-1] for i=1:Ni) for n=1:nbscen)
						);
					end
				end
			end
        end
    end
    return cutviolFlag;
end


###############################################################
###############################################################

#Train model
function train_models_offline_FAD()
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
        xval, thetaval, lb, in_sample = FOSDDP_forward_pass_oneSP_iteration_FAD(lb,xval,thetaval);
        push!(LB,lb);
		#println("LB = ", lb)
        #termination check
        flag, Elapsed = termination_check(iter,relative_gap,LB,start,cutviol_iter);
        if flag != 0
            train_time = Elapsed;
            break;
        end

        #backward pass (if not terminated)
        cutviolFlag = FOSDDP_backward_pass_oneSP_iteration_FAD(lb,xval,thetaval,in_sample);
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
function FOSDDP_eval_offline_FAD()
    start=time();
    OS_paths = Matrix(CSV.read("./data/OOS.csv",DataFrame)); #read the out-of-sample file
	# we do not have this second layer now [REVISION]
	#OS_M = Matrix(CSV.read("./data/inOOS.csv",DataFrame))[:,1] #read the second layer OOS
    objs_fa = zeros(nbOS,Tmin+1);
    
    for s=1:nbOS
        xval = zeros(Ni,Tmin);
        for t=1:Tmin
            #the state is known in the first stage; if not sample a new state k 
            k_t = OS_paths[s,t];
            # we do not have this second layer now [REVISION]
			#m = OS_M[s]; # realization from OS path corresponding to layer 2

            if k_t ∉ absorbing_states
                if t > 1
					#update the RHS
					for i=1:Ni
						set_normalized_rhs(FB1Cons_fa[t,k_t][i], xval[i,t-1]);
						set_normalized_rhs(FB2Cons_fa[t,k_t][i], xval[i,t-1]);
					end
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
                    objs_fa[s,t] = objective_value(m_fa[t,k_t])- value(ϴ_fa[t,k_t]);
                end
            end
        end     
		k_t = OS_paths[s,Tmin+1];
		if k_t in absorbing_states
			continue
		end
		if S[k_t][3] == Nc-1
			# made landfall -> deterministic realization
			for i = 1:Ni
				set_normalized_rhs(FB_final[k_t][i], xval[i,Tmin]);
			end
			for j = 1:Nj
				set_normalized_rhs(dCons_final[k_t][j], SCEN[k_t][j]);
			end
			for i=1:Ni
				for j=1:Nj
					set_objective_coefficient(model_final[k_t], y_final[k_t][i,j], ca[i,j,Tmin]);
				end
			end
			#solve the model
			optimize!(model_final[k_t]);

			#check the status 
			status = termination_status(model_final[k_t]);
			if status != MOI.OPTIMAL
				println("Error in evaluation, final step");
				exit(0);
			else
				#collect values
				objs_fa[s,Tmin+1] = objective_value(model_final[k_t]);
			end
		else
			absorbingT = -1;
			τ = nothing
			for n=1:nbscen
				#identify the period when the hurricane makes landfall 
				τ = findfirst(x -> S[x][3] == Nc-1 && x ∉ absorbing_states, OS_paths[s,:]);
				
				#update the RHS
				if τ === nothing     
					absorbingT = findfirst(x -> S[x][1] == 1, OS_paths[s,:]);
					for i = 1:Ni
						set_normalized_rhs(FB_final[k_t][i], xval[i,Tmin]);
					end
					for j = 1:Nj
						set_normalized_rhs(dCons_final[k_t][j], 0);
					end
					for i=1:Ni
						for j=1:Nj
							set_objective_coefficient(model_final[k_t], y_final[k_t][i,j], ca[i,j,absorbingT]);
						end
					end
					for tt = (Tmin+1):absorbingT 
						objs_fa[s,Tmin+1] += sum(ch[i,tt]*xval[i,Tmin] for i=1:Ni);
					end
				else
					for i = 1:Ni
						set_normalized_rhs(FB_final[k_t][i], xval[i,Tmin]);
					end
					for j = 1:Nj
						set_normalized_rhs(dCons_final[k_t][j], SCEN[OS_paths[s,τ]][j]);
					end
					for i=1:Ni
						for j=1:Nj
							set_objective_coefficient(model_final[k_t], y_final[k_t][i,j], ca[i,j,τ]);
						end
					end
					for tt = (Tmin+1):τ 
						objs_fa[s,Tmin+1] += sum(ch[i,tt]*xval[i,Tmin] for i=1:Ni);
					end
				end
				#solve the subproblem and store the dual information
				optimize!(model_final[k_t]) #solve the model
				status_subproblem = termination_status(model_final[k_t]); #check the status 
				if status_subproblem != MOI.OPTIMAL
					println("Error in Backward Pass, final step");
					println("Model in stage =", t, " and state = ", sample_n, ", in backward pass is ", status);
					exit(0);
				else
					#update the values
					objs_fa[s,Tmin+1] += objective_value(model_final[k_t]);
				end
			end
		end
    end
    fa_bar = mean(sum(objs_fa[:,t] for t=1:(Tmin+1)));
    fa_std = std(sum(objs_fa[:,t] for t=1:(Tmin+1)));
    fa_low = fa_bar-1.96*fa_std/sqrt(nbOS);
    fa_high = fa_bar+1.96*fa_std/sqrt(nbOS);
	println("deterministic FA...");
    println("μ ± 1.96*σ/√NS = ", fa_bar, " ± ", [fa_low,fa_high]);
    elapsed = time() - start;
    return objs_fa, fa_bar, fa_low, fa_high, elapsed
end


