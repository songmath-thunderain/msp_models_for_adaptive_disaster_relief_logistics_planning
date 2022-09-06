s = 1;
t_roll = 1;
x_init = x_0;
#define the model.
master, x, f, θ, subproblem, y2, xCons, dCons, rCons = RH_2SSP_define_models(t_roll,x_init);

#solve the model.
LB_1stRoll, UB_1stRoll, iter_1stRoll, xval_1stRoll, fval_1stRoll, θval_1stRoll = RH_2SSP_solve_roll(s,t_roll,master,subproblem,x,f,θ,y2,xCons,dCons,rCons);

OS_paths = Matrix(CSV.read("./data/OOS.csv",DataFrame)); #read the out-of-sample file
#OS_M = Matrix(CSV.read("./data/inOOS.csv",DataFrame))[:,1] #read the second layer OOS [REVISION: no need anymore]

f1cost = LB_1stRoll-θval_1stRoll;

objs_RH2SSP = zeros(nbOS,T);
objs_RH2SSP[:,1] .= f1cost;

for s=1:nbOS
	#println("s = ", s)
	τ = findfirst(x -> S[x][3] == Nc-1 && x ∉ absorbing_states, OS_paths[s,1:T]);
	x_init = deepcopy(xval_1stRoll[:,1]);
	if τ !== nothing
		# This rolling procedure will stop at t = τ
		for t_roll=2:(τ-1)
			# roll up to t = τ-1
			k_t = OS_paths[s,t_roll]
			if S[k_t][3] == Nc-1 && k_t ∉ absorbing_states
				ξ = SCEN[k_t];
			end

			#define the the model.
			master, x, f, θ, subproblem, y2, xCons, dCons, rCons = RH_2SSP_define_models(t_roll,x_init);
			
			#solve the model.
			LB_Roll, UB_Roll, iter_Roll, xval_Roll, fval_Roll, θval_Roll = RH_2SSP_solve_roll(s,t_roll,master,subproblem,x,f,θ,y2,xCons,dCons,rCons);

			#implement xₜ, pay cost, and pass xₜ to new t+1.
			x_init = deepcopy(xval_Roll[:,1]);
			objs_RH2SSP[s,t_roll] = LB_Roll - θval_Roll;
			if t_roll == (τ-1)
				# Now we get the realization, do the recourse now and finish the rolling procedure
				# TBD
			end
		end
	else
		# This rolling procedure will go all the way to the end
		for t_roll=2:(T-1)
			# roll up to t = T-1
			k_t = OS_paths[s,t_roll]
			if S[k_t][3] == Nc-1 && k_t ∉ absorbing_states
				ξ = SCEN[k_t];
			end

			#define the the model.
			master, x, f, θ, subproblem, y2, xCons, dCons, rCons = RH_2SSP_define_models(t_roll,x_init);
			
			#solve the model.
			LB_Roll, UB_Roll, iter_Roll, xval_Roll, fval_Roll, θval_Roll = RH_2SSP_solve_roll(s,t_roll,master,subproblem,x,f,θ,y2,xCons,dCons,rCons);

			#implement xₜ, pay cost, and pass xₜ to new t+1.
			x_init = deepcopy(xval_Roll[:,1]);
			objs_RH2SSP[s,t_roll] = LB_Roll - θval_Roll;
			if t_roll == (T-1)
				# Now we get the realization, do the recourse now and finish the rolling procedure
				# TBD
			end
		end
	end

RH2SSP_bar = mean(sum(objs_RH2SSP[:,t] for t=1:T));
RH2SSP_std = std(sum(objs_RH2SSP[:,t] for t=1:T));
RH2SSP_low = RH2SSP_bar-1.96*RH2SSP_std/sqrt(nbOS);
RH2SSP_high = RH2SSP_bar+1.96*RH2SSP_std/sqrt(nbOS);
println("μ ± 1.96*σ/√NS = ", RH2SSP_bar, " ± ", [RH2SSP_low,RH2SSP_high]);


fname = "./output/benchmark/roliing2SPresults.csv"
df = CSV.read(fname,DataFrame);

results_RH2SSP = Matrix(df);
results_RH2SSP[inst,1] = 0
results_RH2SSP[inst,2] = RH2SSP_bar
results_RH2SSP[inst,3] = RH2SSP_bar-RH2SSP_low
results_RH2SSP[inst,4] = 0
results_RH2SSP[inst,5] = elapsed_RH2SSP
results_RH2SSP[inst,6] = 0

updf = DataFrame(results_RH2SSP, :auto);
CSV.write(fname,updf)


println("############################################################")
println("############################################################")
