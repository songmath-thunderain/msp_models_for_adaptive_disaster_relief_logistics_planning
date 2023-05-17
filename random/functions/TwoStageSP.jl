# Define the terminal-stage problem: only used when absorbing_option = 1, i.e., MDC/SP operation is allowed to occur
function terminal_model(t_roll,x_init)
 
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
            end
          );

    #######################
    #Define the objective.
    @objective(m,
               Min,
               sum(sum(cb[i,ii,t_roll]*f[i,ii] for ii=1:Ni) for i=1:N0)
              +sum(ch[i,t_roll]*x[i] for i=1:Ni)
              +sum(f[N0,i] for i=1:Ni)*h[t_roll]
              +sum(sum(ca[i,j,t_roll]*y[i,j] for j=1:Nj) for i=1:Ni)
              +sum(z[j] for j=1:Nj)*p
              +sum(v[i] for i=1:Ni)*q       
               );

    #######################
    #Define the constraints.
    dCons = Dict(); #a dictonary to store all the demand constraint
	for i=1:Ni
		@constraint(m, 
					x[i]
				   +sum(f[i,j] for j=1:Ni if j != i)
				   -sum(f[j,i] for j=1:N0 if j != i)
				   +sum(y[i,j] for j=1:Nj)
				   +v[i]
				   == x_init[i]
					);
		@constraint(m, sum(f[i,j] for j=1:Ni if j != i) <= x_init[i]);
	end
	for j=1:Nj
	   dCons[j] = @constraint(m, z[j]+sum(y[i,j] for i=1:Ni) >= 0);
	end

    return m, x, f, y, z, v, dCons
end

#Define first-stage master problem
function RH_2SSP_first_stage(t_roll,k_t,x_init)
    #Note that the static 2SSP corresponds to the case when t_roll = 1 
	nbstages1 = T-t_roll+1;
	if absorbing_option == 0
		# If no MDC/SP operation is allowed in the absorbing state, do not plan for stage T since we know for sure that all states are absorbing
		nbstages1 = T-t_roll;
	end
    m = Model(optimizer_with_attributes(() -> Gurobi.Optimizer(GRB_ENV), "OutputFlag" => 0));
	nbScens = length(nodeScenList[t_roll,k_t]);
	pbScens = nodeScenWeights[t_roll,k_t];
    #######################
    #Define the variables, note that t is a relative index (to the current roll) going from 1 to nbstages1; from relative index to absolute index: t-> t_roll-1+t
    @variables(m,
                begin
                   0 <= x[i=1:Ni,t=1:nbstages1] <= x_cap[i]
                   0 <= f[i=1:N0,ii=1:Ni,t=1:nbstages1]
               -1e8 <= θ[n=1:nbScens] # multicut version! 
                end
              );
    
    #######################
    #Define the objective.
    @objective(m,
               Min,
               sum(sum(sum(cb[i,ii,t_roll-1+t]*f[i,ii,t] for ii=1:Ni) for i=1:N0) for t=1:nbstages1)
              +sum(sum(ch[i,t_roll-1+t]*x[i,t] for i=1:Ni) for t=1:nbstages1)
              +sum(sum(h[t_roll-1+t]*f[N0,i,t] for i=1:Ni) for t=1:nbstages1)
              +sum(θ[n]*pbScens[n] for n=1:nbScens)
               );
    
    #######################
    #Define the constraints.
    for t=1:nbstages1
        for i=1:Ni
            if t == 1
                @constraint(m, x[i,t]
                               +sum(f[i,j,t] for j=1:Ni if j != i)
                               -sum(f[j,i,t] for j=1:N0 if j != i)
                               == x_init[i]
                            ); 
                @constraint(m, sum(f[i,j,t] for j=1:Ni if j != i) <= x_init[i]);                
            else
				@constraint(m, x[i,t-1]
                             -x[i,t]
                             -sum(f[i,j,t] for j=1:Ni if j != i)
                             +sum(f[j,i,t] for j=1:N0 if j != i)==0);
                @constraint(m, sum(f[i,j,t] for j=1:Ni if j != i) <= x[i,t-1]);                
            end
        end
    end
    
    return m, x, f, θ 
    
end


###############################################################
###############################################################

#Define second-stage scenario supbproblem
function RH_2SSP_second_stage()
    #######################
    m = Model(optimizer_with_attributes(() -> Gurobi.Optimizer(GRB_ENV), "OutputFlag" => 0));

    #######################
    #Define the variables.
    @variables(m,
                begin
                   0 <= y[i=1:Ni,j=1:Nj]
                   0 <= z[j=1:Nj]
                   0 <= v[i=1:Ni]
                   0 >= reimbursement
                end
              );
    
    #######################
    #Define the objective.
    @objective(m,
               Min,
              +sum(sum(ca[i,j,1]*y[i,j] for j=1:Nj) for i=1:Ni) # note that we will update the coefficient on y for different scenarios, use ca[i,j,1] as a placeholder for now
              +sum(z[j] for j=1:Nj)*p 
              +sum(v[i] for i=1:Ni)*q 
              +reimbursement
               );

    #######################
    #Define the constraints.
    xCons = Dict(); #a dictonary to store all the inventory constraints
    dCons = Dict(); #a dictonary to store all the demand constraints
    
	for i=1:Ni
		xCons[i] = @constraint(m, sum(y[i,j] for j=1:Nj)+v[i] == 0);             
	end
	for j=1:Nj
	   dCons[j] = @constraint(m, z[j]+sum(y[i,j] for i=1:Ni) == 0);
	end
	rCons = @constraint(m, reimbursement == 0);

    return m, y, xCons, dCons, rCons
end

###############################################################
###############################################################

#defines the two-stage SP models: master problem and subproblem
function RH_2SSP_define_models(t_roll,k_t,x_init)
    #define first stage (master problem) model  
    master, x, f, θ = RH_2SSP_first_stage(t_roll,k_t,x_init);
    
    #define second stage (subproblem) optimality model
    subproblem, y2, xCons, dCons, rCons = RH_2SSP_second_stage();
    
    return master, x, f, θ, subproblem, y2, xCons, dCons, rCons
end


###############################################################
###############################################################


#solves the two-stage SP model
function RH_2SSP_solve_roll(k_t,t_roll,master,subproblem,x,f,θ,y,xCons,dCons,rCons)
# s is the index for the sample path in the out-of-sample test 
# Assumption coming into this function: t_roll must be before the landfall stage of sample path s in the out-of-sample test
	nbstages1 = T-t_roll+1;
	if absorbing_option == 0
		# If no MDC/SP operation is allowed in the absorbing state, do not plan for stage T since we know for sure that all states are absorbing
		nbstages1 = T-t_roll;
	end
	nbScens = length(nodeScenList[t_roll,k_t]);
	pbScens = nodeScenWeights[t_roll,k_t];
	LB = -1e10; 
	UB = 1e10; 
	θval = zeros(nbScens);
    xval = zeros(Ni,nbstages1);
	fval = zeros(N0,Ni,nbstages1);
	solveIter = 0;
	flag = 1;
#while flag == 1
	while (UB-LB)*1.0/max(1e-10,abs(LB)) > ϵ && abs(UB-LB) > ϵ && solveIter < 100 # WARNING: Temporary fix here!
        # solve first stage
		flag = 0;
		solveIter = solveIter + 1;
        LB, xval, fval, θval = solve_first_stage(LB,xval,fval,θval,master,x,f,θ);
		firstCost = LB - sum(θval[n]*pbScens[n] for n=1:nbScens);
		# solve second stage 
		flag, Qbar = solve_second_stage(k_t,t_roll,xval,fval,θval,master,subproblem,x,f,θ,y,xCons,dCons,rCons);
		if flag != -1
			UB = min(firstCost+Qbar,UB);
		end
    end
	if solveIter == 100
		println("# iterations is maxed out!");
		print("LB = ", LB);
		println(", UB = ", UB);
	end
    return LB, UB, xval, fval, θval
end


###############################################################
###############################################################


#solves the first-stage problem
function solve_first_stage(LB,xval,fval,θval,master,x,f,θ)
    optimize!(master); #solve the model
    status_master = termination_status(master); #check the status 
    if status_master != MOI.OPTIMAL
        println("Master problem status is: ", status_master, " === Oops! :/");
        return
        #exit(0);
    else
        #update the values
        LB = objective_value(master);
        xval = value.(x);
        fval = value.(f);
        θval = value.(θ);      
    end

	return LB, xval, fval, θval
end

###############################################################
###############################################################

#solves the second-stage problem
function solve_second_stage(k_t,t_roll,xval,fval,θval,master,subproblem,x,f,θ,y,xCons,dCons,rCons)
	nbScens = length(nodeScenList[t_roll,k_t]);
	pbScens = nodeScenWeights[t_roll,k_t];
    flag = 0;
    nbstages1 = T-t_roll+1;
	if absorbing_option == 0
		# If no MDC/SP operation is allowed in the absorbing state, do not plan for stage T since we know for sure that all states are absorbing
		nbstages1 = T-t_roll;
	end
    Q = zeros(nbScens); #list for all the optimal values
    pi1 = Array{Any,1}(undef,nbScens); #list for all the dual multiplies of the first set of constraints
    pi2 = Array{Any,1}(undef,nbScens); #list for all the dual multiplies of the second set of constraints
    pi3 = zeros(nbScens); # dual multipliers of the third set of constraints
    Qbar = 0;

    for n=1:nbScens
        #identify the period when the hurricane makes landfall 
		RH_2SSP_update_RHS(nodeScenList[t_roll,k_t][n][1],nodeScenList[t_roll,k_t][n][2],subproblem,xCons,dCons,rCons,xval,fval,y,t_roll);
        #solve the subproblem and store the dual information
        Q[n], pi1[n], pi2[n], pi3[n], flag = solve_scen_subproblem(subproblem,xCons,dCons,rCons);
		if flag == -1
            println("subproblem status is infeasible?!");
            exit(0);
		end
    end
    Qbar = sum(Q[n]*pbScens[n] for n=1:nbScens);

    # cut generation: multi-cut version
    for n=1:nbScens
	    if (Q[n]-θval[n])/max(1e-10,abs(Q[n])) > ϵ && Q[n]-θval[n] > ϵ
			#println("cutviol[", n, "] = ", Q[n]-θval[n], "Q[", n, "] = ", Q[n], "θval[", n, "] = ", θval[n]);
			tt = nodeScenList[t_roll,k_t][n][1];
			# tt is the terminal stage
			if absorbing_option == 0
				@constraint(master,
				θ[n]-sum(pi1[n][i]*x[i,tt-t_roll] for i=1:Ni)-pi3[n]*(-sum((sum(sum(cb[i,ii,t_roll+t-1]*f[i,ii,t] for ii=1:Ni) for i=1:N0)+sum(ch[i,t_roll+t-1]*x[i,t] for i=1:Ni)+sum(f[N0,i,t] for i=1:Ni)*h[t_roll+t-1]) for t = (tt+1-t_roll):nbstages1)) 
					>= 
				 Q[n]-sum(pi1[n][i]*xval[i,tt-t_roll] for i=1:Ni)-pi3[n]*(-sum((sum(sum(cb[i,ii,t_roll+t-1]*fval[i,ii,t] for ii=1:Ni) for i=1:N0)+sum(ch[i,t_roll+t-1]*xval[i,t] for i=1:Ni)+sum(fval[N0,i,t] for i=1:Ni)*h[t_roll+t-1]) for t = (tt+1-t_roll):nbstages1))
					);
			else
				@constraint(master,
				θ[n]-sum(pi1[n][i]*x[i,tt-t_roll+1] for i=1:Ni)-pi3[n]*(-sum((sum(sum(cb[i,ii,t_roll+t-1]*f[i,ii,t] for ii=1:Ni) for i=1:N0)+sum(ch[i,t_roll+t-1]*x[i,t] for i=1:Ni)+sum(f[N0,i,t] for i=1:Ni)*h[t_roll+t-1]) for t = (tt+2-t_roll):nbstages1)) 
					>= 
				 Q[n]-sum(pi1[n][i]*xval[i,tt-t_roll+1] for i=1:Ni)-pi3[n]*(-sum((sum(sum(cb[i,ii,t_roll+t-1]*fval[i,ii,t] for ii=1:Ni) for i=1:N0)+sum(ch[i,t_roll+t-1]*xval[i,t] for i=1:Ni)+sum(fval[N0,i,t] for i=1:Ni)*h[t_roll+t-1]) for t = (tt+2-t_roll):nbstages1))
					);
			end
			flag = 1;
	    end
    end
    return flag, Qbar
end


###############################################################
###############################################################
#solves scenario subproblem of the second stage
function solve_scen_subproblem(subproblem,xCons,dCons,rCons)
    flag = 0
    optimize!(subproblem) #solve the model
    status_subproblem = termination_status(subproblem); #check the status 
    Qtemp = 0;
    pi1temp = zeros(Ni);
    pi2temp = zeros(Nj);
    pi3temp = 0;
    if status_subproblem != MOI.OPTIMAL
        if status_subproblem == MOI.INFEASIBLE
            flag = -1;
        end
    else
        #update the values
        Qtemp = objective_value(subproblem);
		pi1temp = zeros(Ni);
        pi2temp = zeros(Nj);
        for i=1:Ni
            pi1temp[i] = shadow_price(xCons[i]);
        end
		for j=1:Nj
            pi2temp[j] = shadow_price(dCons[j]);            
        end
		pi3temp = shadow_price(rCons);
    end
    return Qtemp, pi1temp, pi2temp, pi3temp, flag
end

###############################################################
###############################################################

#updates the RHS of the 2nd-stage constraints and objective coefficients
function RH_2SSP_update_RHS(absorbingT,k_t,subproblem,xCons,dCons,rCons,xval,fval,y,t_roll)
	nbstages1 = T-t_roll+1;
	if absorbing_option == 0
		# If no MDC/SP operation is allowed in the absorbing state, do not plan for stage T since we know for sure that all states are absorbing
		nbstages1 = T-t_roll;
	end

	for i=1:Ni
		if absorbing_option == 0
			set_normalized_rhs(xCons[i],xval[i,absorbingT-t_roll]);	
		else
			set_normalized_rhs(xCons[i],xval[i,absorbingT-t_roll+1]);
		end
	end

	for j=1:Nj
	    if S[k_t][1] != 1           
        	set_normalized_rhs(dCons[j], SCEN[k_t][j]);
        else
            set_normalized_rhs(dCons[j], 0);
        end
	end

	if absorbingT == T
		# Plan exactly until the landfall time -- no reimbursement occurred!
		set_normalized_rhs(rCons,0);
	else
		updatedRHS = 0;
		if absorbing_option == 0
			# reimburse the operational cost starting from the terminal stage, since the terminal stage does not allow operation
			updatedRHS = -sum((sum(sum(cb[i,ii,t_roll+t-1]*fval[i,ii,t] for ii=1:Ni) for i=1:N0)+sum(ch[i,t_roll+t-1]*xval[i,t] for i=1:Ni)+sum(fval[N0,i,t] for i=1:Ni)*h[t_roll+t-1]) for t = (absorbingT+1-t_roll):nbstages1);
		else
			# reimburse the operational cost if they occur after the terminal stage: starting from stage (τ+1)-t_roll+1
			updatedRHS = -sum((sum(sum(cb[i,ii,t_roll+t-1]*fval[i,ii,t] for ii=1:Ni) for i=1:N0)+sum(ch[i,t_roll+t-1]*xval[i,t] for i=1:Ni)+sum(fval[N0,i,t] for i=1:Ni)*h[t_roll+t-1]) for t = (absorbingT+2-t_roll):nbstages1); 
		end
		set_normalized_rhs(rCons,updatedRHS);
	end
	# Also need to update the coefficients of y[i,j] variables in the 2nd stage
	for i=1:Ni
		for j=1:Nj
			set_objective_coefficient(subproblem, y[i,j], ca[i,j,absorbingT]);
		end
	end
end

###############################################################
###############################################################
