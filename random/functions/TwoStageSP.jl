#Define first-stage master problem
function RH_2SSP_first_stage(t_roll,x_init)
    #Note that the static 2SSP corresponds to the case when t_roll = 1 
	nbstages1 = T-t_roll+1;
    m = Model(optimizer_with_attributes(() -> Gurobi.Optimizer(GRB_ENV), "OutputFlag" => 0));

    #######################
    #Define the variables, note that t is a relative index (to the current roll) going from 1 to nbstages1; from relative index to absolute index: t-> t_roll-1+t
    @variables(m,
                begin
                   0 <= x[i=1:Ni,t=1:nbstages1] <= x_cap[i]
                   0 <= f[i=1:N0,ii=1:Ni,t=1:nbstages1] <= f_cap[i,ii]
               -1e8 <= θ[n=1:nbscen] # multicut version! 
                end
              );
    
    #######################
    #Define the objective.
    @objective(m,
               Min,
               sum(sum(sum(cb[i,ii,t_roll-1+t]*f[i,ii,t] for ii=1:Ni) for i=1:N0) for t=1:nbstages1)
              +sum(sum(ch[i,t_roll-1+t]*x[i,t] for i=1:Ni) for t=1:nbstages1)
              +sum(sum(h[t_roll-1+t]*f[N0,i,t] for i=1:Ni) for t=1:nbstages1)
              +sum(θ[n]*1.0/nbscen for n=1:nbscen)
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
    m = Model(optimizer_with_attributes(() -> Gurobi.Optimizer(GRB_ENV), "OutputFlag" => 0, "Presolve" => 0));

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
function RH_2SSP_define_models(t_roll,x_init)
    #define first stage (master problem) model  
    master, x, f, θ = RH_2SSP_first_stage(t_roll,x_init);
    
    #define second stage (subproblem) optimality model
    subproblem, y2, xCons, dCons, rCons = RH_2SSP_second_stage();
    
    return master, x, f, θ, subproblem, y2, xCons, dCons, rCons
end


###############################################################
###############################################################

#initialize parameter for two-stage model
function initialize(s,t_roll)
# s is the index for the sample path in the out-of-sample test 
    nbstages1 = T-t_roll+1;
	LB = -1e10; 
	UB = 1e10; 
	θval = zeros(nbscen);
    xval = zeros(Ni,nbstages1);
	fval = zeros(N0,Ni,nbstages1);
	# Do sampling to create (in-sample) scenarios 
	osfname = "./data/OOS"*string(k_init)*".csv";
    allscen = Matrix(CSV.read(osfname,DataFrame)); #read the out-of-sample file
	# In the first roll, always choose the nbscen scenarios out of the total of 10000 OOS once every 10000/nbscen 
    scen = allscen[collect(1:convert(Int,10000/nbscen):10000),1:T]; # note that this is just an initialization for scen

	# when t_roll = 1, i.e., the first roll, it will always be the evenly selected scenarios from allscen

	# In later rolls, create in-sample scenarios by sampling
    if t_roll > 1
        for n=1:nbscen
            for t=2:t_roll
                scen[n,t] = allscen[s,t];
            end
            
            for t=(t_roll+1):T
                scen[n,t] = MC_sample(scen[n,t-1]);
            end
        end
    end
    qprob = fill(1/nbscen,nbscen);
	return LB, UB, xval, fval, θval, scen, qprob
end

###############################################################
###############################################################

#initialize parameter for two-stage model
function initialize2(kk,t_roll)
# kk is the starting state for the two-stage SP model: this is to accommodate the wait-and-see heuristic approach    
	nbstages1 = T-t_roll+1;
	LB = -1e10; 
	UB = 1e10; 
	θval = zeros(nbscen);
    xval = zeros(Ni,nbstages1);
	fval = zeros(N0,Ni,nbstages1);
	# Do sampling to create (in-sample) scenarios 
	osfname = "./data/OOS"*string(k_init)*".csv";
    allscen = Matrix(CSV.read(osfname,DataFrame)); #read the out-of-sample file
	# In the first roll, always choose the nbscen scenarios out of the total of 10000 OOS once every 10000/nbscen 
    scen = allscen[collect(1:convert(Int,10000/nbscen):10000),1:T]; # note that this is just an initialization for scen

	# when t_roll = 1, i.e., the first roll, it will always be the evenly selected scenarios from allscen

	# In later rolls, create in-sample scenarios by sampling
    if t_roll > 1
        for n=1:nbscen
            for t=2:t_roll
                scen[n,t] = allscen[1,1]; # prior to t_roll, just choose arbitrarily 
            end
            scen[n,t_roll] = kk;
            for t=(t_roll+1):T
                scen[n,t] = MC_sample(scen[n,t-1]);
            end
        end
    end
    qprob = fill(1/nbscen,nbscen);
	return LB, UB, xval, fval, θval, scen, qprob
end

###############################################################
###############################################################


#solves the two-stage SP model
function RH_2SSP_solve_roll(s,t_roll,master,subproblem,x,f,θ,y,xCons,dCons,rCons)
# s is the index for the sample path in the out-of-sample test 
# Assumption coming into this function: t_roll must be before the landfall stage of sample path s in the out-of-sample test
    nbstages1 = T-t_roll+1;
    LB, UB, xval, fval, θval, scen, qprob = initialize(s,t_roll);
	solveIter = 0;
    while (UB-LB)*1.0/max(1e-10,abs(LB)) > ϵ && abs(UB-LB) > ϵ && solveIter < 100 # WARNING: Temporary fix here!
        # solve first stage
		solveIter = solveIter + 1;
        LB, xval, fval, θval = solve_first_stage(LB,xval,fval,θval,master,x,f,θ);
		firstCost = LB - sum(θval[n]*qprob[n] for n=1:nbscen);
        if t_roll < T
            # solve second stage 
            flag, Qbar = solve_second_stage(t_roll,xval,fval,θval,scen,qprob,master,subproblem,x,f,θ,y,xCons,dCons,rCons);
            if flag != -1
                UB = min(firstCost+Qbar,UB);
            end
        else
            println("exiting here");
            break;
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
###############################################################\

#solves the two-stage SP model
function RH_2SSP_solve_roll2(kk,t_roll,master,subproblem,x,f,θ,y,xCons,dCons,rCons)
# kk is the starting state for the two-stage SP model: this is to accommodate the wait-and-see heuristic approach 
# Assumption coming into this function: t_roll must be before the landfall stage of sample path s in the out-of-sample test
    nbstages1 = T-t_roll+1;
    LB, UB, xval, fval, θval, scen, qprob = initialize2(kk,t_roll);
	solveIter = 0;
    while (UB-LB)*1.0/max(1e-10,abs(LB)) > ϵ && abs(UB-LB) > ϵ && solveIter < 100 # WARNING: Temporary fix here!
        # solve first stage
		solveIter = solveIter + 1;
        LB, xval, fval, θval = solve_first_stage(LB,xval,fval,θval,master,x,f,θ);
		firstCost = LB - sum(θval[n]*qprob[n] for n=1:nbscen);
        if t_roll < T
            # solve second stage 
            flag, Qbar = solve_second_stage(t_roll,xval,fval,θval,scen,qprob,master,subproblem,x,f,θ,y,xCons,dCons,rCons);
            if flag != -1
                UB = min(firstCost+Qbar,UB);
            end
        else
            println("exiting here");
            break;
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
function solve_second_stage(t_roll,xval,fval,θval,scen,qprob,master,subproblem,x,f,θ,y,xCons,dCons,rCons)
    flag = 0;
    nbstages1 = T-t_roll+1;
    Q = zeros(nbscen); #list for all the optimal values
    pi1 = Array{Any,1}(undef,nbscen); #list for all the dual multiplies of the first set of constraints
    pi2 = Array{Any,1}(undef,nbscen); #list for all the dual multiplies of the second set of constraints
    pi3 = zeros(nbscen); # dual multipliers of the third set of constraints
    Qbar = 0;
    absorbingT = -1;
    τ = nothing

    for n=1:nbscen
        #identify the period when the hurricane makes landfall 
        τ = findfirst(x -> (S[x][3] == Nc && S[x][1] != 1), scen[n,:]);
        
        #update the RHS
        if τ === nothing     
	    	absorbingT = findfirst(x -> S[x][1] == 1, scen[n,:]);
			#print("n = ", n);
			#print(", τ = ", τ);
			#print(", absorbingT = ", absorbingT);
			#println(", t_roll = ", t_roll);
			#if absorbingT < t_roll
			#	absorbingT = t_roll; # WARNING: temporary fix!
			#end
            RH_2SSP_update_RHS(absorbingT,scen[n,absorbingT],subproblem,xCons,dCons,rCons,xval,fval,y,t_roll);
        else
			#print("n = ", n);
			#print(", τ = ", τ);
			#println(", t_roll = ", t_roll);
	    	absorbingT = -1;
			#if τ < t_roll
			#	τ = t_roll; # WARNING: temporary fix!
			#end
            RH_2SSP_update_RHS(τ,scen[n,τ],subproblem,xCons,dCons,rCons,xval,fval,y,t_roll);
        end
        
        #solve the subproblem and store the dual information
        Q[n], pi1[n], pi2[n], pi3[n], flag = solve_scen_subproblem(subproblem,xCons,dCons,rCons);

		if flag == -1
            println("subproblem status is infeasible?!");
            exit(0);
		end
    end
    Qbar = sum(Q[n]*qprob[n] for n=1:nbscen);

    # cut generation: multi-cut version
    for n=1:nbscen
	    if (Q[n]-θval[n])/max(1e-10,abs(Q[n])) > ϵ && abs(Q[n]-θval[n]) > ϵ
			τ = findfirst(x -> (S[x][3] == Nc && S[x][1] != 1), scen[n,:]);
			tt = -1;
			if τ === nothing
		    	tt = findfirst(x -> S[x][1] == 1, scen[n,:]);
			else
		    	tt = τ;
			end
			if tt < t_roll
				tt = t_roll; # WARNING: temporary fix!
			end	
			@constraint(master,
			θ[n]-sum(pi1[n][i]*x[i,tt-t_roll+1] for i=1:Ni)-pi3[n]*(-sum((sum(sum(cb[i,ii,t_roll+t-1]*f[i,ii,t] for ii=1:Ni) for i=1:N0)+sum(ch[i,t_roll+t-1]*x[i,t] for i=1:Ni)+sum(f[N0,i,t] for i=1:Ni)*h[t_roll+t-1]) for t = (tt+2-t_roll):nbstages1)) 
				>= 
			 Q[n]-sum(pi1[n][i]*xval[i,tt-t_roll+1] for i=1:Ni)-pi3[n]*(-sum((sum(sum(cb[i,ii,t_roll+t-1]*fval[i,ii,t] for ii=1:Ni) for i=1:N0)+sum(ch[i,t_roll+t-1]*xval[i,t] for i=1:Ni)+sum(fval[N0,i,t] for i=1:Ni)*h[t_roll+t-1]) for t = (tt+2-t_roll):nbstages1))
				);
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
function RH_2SSP_update_RHS(τ,k_t,subproblem,xCons,dCons,rCons,xval,fval,y,t_roll)
	# Note that the assumption coming into this function is that τ > t_roll if τ is not nothing
	nbstages1 = T-t_roll+1;
	for i=1:Ni
	    if τ === nothing
			println("We shouldn't come here any more!");
			exit(0);
			#set_normalized_rhs(xCons[i],xval[i,T-t_roll+1]);
	    else
			set_normalized_rhs(xCons[i],xval[i,τ-t_roll+1]);
	    end
	end

	for j=1:Nj
	    if S[k_t][1] != 1           
        	set_normalized_rhs(dCons[j], SCEN[k_t][j]);
        else
            set_normalized_rhs(dCons[j], 0);
        end
	end

	if τ === nothing 
	    # nothing to reimburse here
        println("We shouldn't come here any more!");
	    exit(0);
	    #set_normalized_rhs(rCons, 0)
    else
		if τ == T
			# Plan exactly until the landfall time -- no reimbursement occurred!
			set_normalized_rhs(rCons,0);
		else
			updatedRHS = -sum((sum(sum(cb[i,ii,t_roll+t-1]*fval[i,ii,t] for ii=1:Ni) for i=1:N0)
						+sum(ch[i,t_roll+t-1]*xval[i,t] for i=1:Ni)  
						+sum(fval[N0,i,t] for i=1:Ni)*h[t_roll+t-1]) for t = (τ+2-t_roll):nbstages1); 
			set_normalized_rhs(rCons,updatedRHS);
		end
	    # Also need to update the coefficients of y[i,j] variables in the 2nd stage
	    for i=1:Ni
	        for j=1:Nj
			    set_objective_coefficient(subproblem, y[i,j], ca[i,j,τ]);
			end
	    end
    end
end

###############################################################
###############################################################
