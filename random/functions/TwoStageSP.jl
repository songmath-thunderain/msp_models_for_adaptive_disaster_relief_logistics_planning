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
               -1e8 <= θ 
                end
              );
    
    #######################
    #Define the objective.
    @objective(m,
               Min,
               sum(sum(sum(cb[i,ii,t_roll-1+t]*f[i,ii,t] for ii=1:Ni) for i=1:N0) for t=1:nbstages1)
              +sum(sum(ch[i,t_roll-1+t]*x[i,t] for i=1:Ni) for t=1:nbstages1)
              +sum(sum(f[N0,i,t] for i=1:Ni)*h[t_roll-1+t] for t=1:nbstages1)
              +θ 
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
              +sum(sum(ca[i,j,1]*y[i,j] for j=1:Nj) for i=1:Ni) # note that we will update the coefficient on y for different scenarios
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
	rCons = @constraint(m, reimbursement == 0)

    return m, y, xCons, dCons, rCons
end

###############################################################
###############################################################
#=
# No need to define a feasibility problem in this revised version.
function feasability_problem(t_roll,nbstages1,nbstages2)
    #######################
    #Define the model.
    m = Model(optimizer_with_attributes(() -> Gurobi.Optimizer(GRB_ENV), "OutputFlag" => 0));

    #######################
    #Define the variables.
    @variables(m,
                begin
                   0 <= y[i=1:Ni,j=1:Nj,t=1:nbstages2]
                   0 <= z[j=1:Nj,t=1:nbstages2]
                   0 <= v[i=1:Ni,t=1:nbstages2]
                   0 >= reimbursement[t=1:nbstages2]
                   0 <= a_plus[i=1:Ni,t=1:nbstages2]
                   0 <= a_minus[i=1:Ni,t=1:nbstages2]
                end
              );
    
    #######################
    #Define the objective.
    @objective(m,
               Min,
               sum(sum(a_plus[i,t] for i=1:Ni) for t=1:nbstages2)
              +sum(sum(a_minus[i,t] for i=1:Ni) for t=1:nbstages2)
               );

    #######################
    #Define the constraints.
    ConsFB = Dict(); #a dictonary to store first set of constraints
    dCons = Dict(); #a dictonary to store all the demand constraint
    rCons = Dict(); #a dictonary to store all the reimbursement constraint
    
    for t=1:nbstages2
        for i=1:Ni
            ConsFB[t,i] = @constraint(m, sum(y[i,j,t] for j=1:Nj)+v[i,t]+a_plus[i,t]-a_minus[i,t]== 0);         
        end
        for j=1:Nj
           dCons[t,j] = @constraint(m, z[j,t]+sum(y[i,j,t] for i=1:Ni) >= 0);
        end
        rCons[t] = @constraint(m, reimbursement[t] == 0)
    end

    return m, y, z, v, ConsFB, dCons, rCons
end
=#
###############################################################
###############################################################

#defines the two-stage SP models: master problem and subproblem
function RH_2SSP_define_models(t_roll,x_init)
    #define first stage (master problem) model  
    master, x, f, θ = RH_2SSP_first_stage(t_roll,x_init)
    
    #define second stage (subproblem) optimality model
    subproblem, y2, xCons, dCons, rCons = RH_2SSP_second_stage()
    
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
	θval = 0;
    xval = zeros(Ni,nbstages1);
	fval = zeros(N0,Ni,nbstages1);
	# Do sampling to create (in-sample) scenarios 
    allscen = Matrix(CSV.read("./data/OOS.csv",DataFrame));
    scen = allscen[collect(1:convert(Int,10000/nbscen):10000),1:T]; # note that this is just an initialization for scen

    if t_roll > 1
        for n=1:nbscen
            for t=2:t_roll
                scen[n,t] = allscen[s,t]
            end
            
            for t=t_roll+1:T
                scen[n,t] = MC_sample(scen[n,t-1]);
            end
        end
    end
    qprob = fill(1/nbscen,nbscen)
	return LB, UB, xval, fval, θval, scen, qprob
end

###############################################################
###############################################################

#solves the two-stage SP model
function RH_2SSP_solve_roll(s,t_roll,master,subproblem,x,f,θ,y,xCons,dCons,rCons)
# s is the index for the sample path in the out-of-sample test 
# Assumption coming into this function: t_roll must be before the landfall stage of sample path s in the out-of-sample test
    nbstages1 = T-t_roll+1;
    LB, UB, xval, fval, θval, scen, qprob = initialize(s,t_roll)
    while (UB-LB)*1.0/max(1e-10,abs(LB)) > ϵ 
        # solve first stage
        solve_first_stage(LB,xval,fval,θval,master,x,f,θ);
        if t_roll < T
            # solve second stage 
            flag, Qbar = solve_second_stage(t_roll,xval,fval,θval,scen,qprob,master,subproblem,x,f,θ,y,xCons,dCons,rCons)
            if flag != -1
                UB = min(LB-θval+Qbar,UB);
            end
        else
            println("exiting here")
            break
        end
    end

    return LB, UB, xval, fval, θval
end


###############################################################
###############################################################

#solves the first-stage problem
function solve_first_stage(LB,xval,fval,θval,master,x,f,θ)
    optimize!(master) #solve the model
    status_master = termination_status(master); #check the status 
    if status_master != MOI.OPTIMAL
        println("Master problem status is: ", status_master, " === Oops! :/")
        return
        #exit(0);
    else
        #update the values
        LB = objective_value(master);
        xval = value.(x);
        fval = value.(f);
        θval = value(θ);      
    end
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
    
    τ = nothing
    for n=1:nbscen
        #identify the period when the hurricane makes landfall 
        τ = findfirst(x -> S[x][3] == Nc-1 && x ∉ absorbing_states, scen[n,:]);
        
        #update the RHS
        if τ === nothing     
            RH_2SSP_update_RHS(τ,scen[n,end],subproblem,xCons,dCons,rCons,xval,fval,y,t_roll)
        else
            RH_2SSP_update_RHS(τ,scen[n,τ],subproblem,xCons,dCons,rCons,xval,fval,y,t_roll)
        end
        
        #solve the subproblem and store the dual information
        flag = solve_scen_subproblem(Q,pi1,pi2,pi3,n,subproblem,xCons,dCons,rCons)
        if flag == -1
        	println("subproblem status is infeasible?!")
        	exit(0);
		end
    end
	Qbar = sum(Q[n]*qprob[n] for n=1:nbscen);
	#add a cut (if any)
	
	if (Qbar-θval)/max(1e-10,abs(Qbar)) > ϵ
		if τ === nothing
			# no reimbursement
			@constraint(master,
			θ-sum(qprob[n]*sum(pi1[n][i]*x[i,nbstages1] for i=1:Ni) for n=1:nbscen) 
			>= 
		 Qbar-sum(qprob[n]*sum(pi1[n][i]*xval[i,nbstages1] for i=1:Ni) for n=1:nbscen)
			);
		else
			# has reimbursement
			@constraint(master,
			θ-sum(qprob[n]*sum(pi1[n][i]*x[i,τ-t_roll+1] for i=1:Ni) for n=1:nbscen)-sum(qprob[n]*pi3[n]*(-sum(sum(sum(cb[i,ii,t_roll+t-1]*f[i,ii,t] for ii=1:Ni) for i=1:N0)+sum(ch[i,t_roll+t-1]*x[i,t] for i=1:Ni)+sum(f[N0,i,t] for i=1:Ni)*h[t_roll+t-1] for t = (τ+2-t_roll):nbstages1)) for n=1:nbscen) 
			>= 
		 Qbar-sum(qprob[n]*sum(pi1[n][i]*xval[i,τ-t_roll+1] for i=1:Ni) for n=1:nbscen)-sum(qprob[n]*pi3[n]*(-sum(sum(sum(cb[i,ii,t_roll+t-1]*fval[i,ii,t] for ii=1:Ni) for i=1:N0)+sum(ch[i,t_roll+t-1]*xval[i,t] for i=1:Ni)+sum(fval[N0,i,t] for i=1:Ni)*h[t_roll+t-1] for t = (τ+2-t_roll):nbstages1)) for n=1:nbscen)
			);
		end
		flag = 1;
	end
    return flag, Qbar
end


###############################################################
###############################################################
#solves scenario subproblem of the second stage
function solve_scen_subproblem(Q,pi1,pi2,pi3,n,subproblem,xCons,dCons,rCons)
    flag = 0
    optimize!(subproblem) #solve the model
    status_subproblem = termination_status(subproblem); #check the status 
    if status_subproblem != MOI.OPTIMAL
        if status_subproblem == MOI.INFEASIBLE
            flag = -1
        end
    else
        #update the values
        Q[n] = objective_value(subproblem);
		pi1temp = zeros(Ni);
        pi2temp = zeros(Nj);
        for i=1:Ni
            pi1temp[i] = shadow_price(xCons[i]);
        end
		for j=1:Nj
            pi2temp[j] = shadow_price(dCons[j]);            
        end
		pi1[n] = pi1temp;
		pi2[n] = pi2temp;
		pi3[n] = shadow_price(rCons);
    end
    return flag
end

###############################################################
###############################################################

#updates the RHS of the 2nd-stage constraints and objective coefficients
function RH_2SSP_update_RHS(τ,k_t,subproblem,xCons,dCons,rCons,xval,fval,y,t_roll)
	# Note that the assumption coming into this function is that τ > t_roll if τ is not nothing
	nbstages1 = T-t_roll+1;
	for i=1:Ni
		if τ === nothing
			set_normalized_rhs(xCons[i],xval[i,T-t_roll+1]);
		else
			set_normalized_rhs(xCons[i],xval[i,τ-t_roll+1]);
		end
	end

	for j=1:Nj
		if k_t ∉ absorbing_states           
            set_normalized_rhs(dCons[j], SCEN[k_t][j]);
        else
            set_normalized_rhs(dCons[j], 0);
        end
	end

	if τ === nothing 
		# nothing to reimburse here
        set_normalized_rhs(rCons, 0)
    else
        set_normalized_rhs(rCons,
                    -sum(sum(sum(cb[i,ii,t_roll+t-1]*fval[i,ii,t] for ii=1:Ni) for i=1:N0)
                    +sum(ch[i,t_roll+t-1]*xval[i,t] for i=1:Ni)  
                    +sum(fval[N0,i,t] for i=1:Ni)*h[t_roll+t-1] for t = (τ+2-t_roll):nbstages1) 
                    );
    end

	# Also need to update the coefficients of y[i,j] variables in the 2nd stage
	if τ != nothing
		for i=1:Ni
			for j=1:Nj
				set_objective_coefficient(subproblem, y[i,j], ca[i,j,τ]);
			end
		end
	end
end

###############################################################
###############################################################
#=
# No need to generate a feasability cut in this revised version.
function generate_feasability_cut(feas_m,feas_y,feas_z,feas_v,feas_ConsFB,feas_dCons,feas_rCons,n,scen,nbstages1,xval,fval,t_roll,nbstages2,master,x,f)
    #identify the period where the hurricane makes landfall 
    τ = findfirst(x -> S[x][3] == Nc-1 && x ∉ absorbing_states, scen[n,:]);

    #update the RHS
    if τ === nothing
        RH_2SSP_update_RHS(τ,scen[n,end],nbstages1,feas_ConsFB,feas_dCons,feas_rCons,xval,fval,t_roll)
    else
        RH_2SSP_update_RHS(τ,scen[n,τ],nbstages1,feas_ConsFB,feas_dCons,feas_rCons,xval,fval,t_roll)
    end

    optimize!(feas_m) #solve the model
    status_feasability = termination_status(feas_m); #check the status 
    if status_feasability != MOI.OPTIMAL
        println("Feasability problem status is: ", status_feasability, " === Oops! :/")
        println("subproblem = ", feas_m)
        exit(0);
    else
        #update the values
        Q = objective_value(feas_m);
        pi1 = zeros(nbstages2,Ni)
        pi2 = zeros(nbstages2)
        
        for t=1:nbstages2
            for i=1:Ni
                pi1[t,i] = dual(feas_ConsFB[t,i]);
            end
            pi2[t] = dual(feas_rCons[t]);            
        end
        if τ === nothing
            @constraint(master,
            0>=sum(pi1[t-1,i]*(x[i,t-1]-x[i,t]-sum(f[i,j,t] for j=1:Ni if j != i)+sum(f[j,i,t] for j=1:N0 if j != i)) for i=1:Ni, t=2:nbstages1)
            );
        else
            @constraint(master,        
            0>=sum(pi1[t-1,i]*(x[i,t-1]-x[i,t]-sum(f[i,j,t] for j=1:Ni if j != i)+sum(f[j,i,t] for j=1:N0 if j != i)) for i=1:Ni, t=2:nbstages1)+sum(pi2[t-1]*(sum(sum(cb[i,ii,t_roll+t-1]*f[i,ii,t] for ii=1:Ni) for i=1:N0)+sum(ch[i,t_roll+t-1]*x[i,t] for i=1:Ni)+sum(f[N0,i,t] for i=1:Ni)*h[t_roll+t-1]) for t=2:nbstages1 if t>τ));
        end
    end
end
=#
