s = 1;
t_roll = 1;
x_init = x_0;

start = time();

#define the model.
master, x, f, θ, subproblem, y2, xCons, dCons, rCons = RH_2SSP_define_models(t_roll,x_init);

#solve the model.
LB_1stRoll, UB_1stRoll, xval_1stRoll, fval_1stRoll, θval_1stRoll = RH_2SSP_solve_roll(s,t_roll,master,subproblem,x,f,θ,y2,xCons,dCons,rCons);
objGo = UB_1stRoll; # objGo: expected objval if implementing a two-stage plan now
Go = false;

# Do sampling to create (in-sample) scenarios 
osfname = "./data/OOS"*string(k_init)*".csv";
allscen = Matrix(CSV.read(osfname,DataFrame)); #read the out-of-sample file

let objNoGo = 0.0; # objNoGo: expected obj if wait to implement a two-stage plan in the next stage
	curr_state = allscen[s,t_roll];
	# sample all the next-stage states and see the corresponding expected objval if implementing a two-stage plan then 
	for kk = 1:K
		if P_joint[curr_state, kk] > 1e-8
			# plan with starting state of kk, initial condition of x_init, and t_roll
			master_nextState, x_nextState, f_nextState, θ_nextState, subproblem_nextState, y2_nextState, xCons_nextState, dCons_nextState, rCons_nextState = RH_2SSP_define_models(t_roll+1,x_init);
			LB_nextState, UB_nextState, xval_nextState, fval_nextState, θval_nextState = RH_2SSP_solve_roll2(kk,t_roll+1,master_nextState,subproblem_nextState,x_nextState,f_nextState,θ_nextState,y2_nextState,xCons_nextState,dCons_nextState,rCons_nextState);
			objNoGo += P_joint[curr_state,kk]*UB_nextState;
		end
	end
	if objGo < objNoGo
		Go = true;
	end
end

if Go == true
	println("same as the static two-stage!");
	exit(0);
end
	
OS_paths = Matrix(CSV.read(osfname,DataFrame)); #read the out-of-sample file
#OS_M = Matrix(CSV.read("./data/inOOS.csv",DataFrame))[:,1] #read the second layer OOS [REVISION: no need anymore]

objs_OOS = zeros(nbOS);

for s=1:nbOS
	τ = findfirst(x -> S[x][3] == Nc-1 && x ∉ absorbing_states, OS_paths[s,1:T]);
	if τ === nothing
		absorbingT = findfirst(x -> S[x][1] == 1, OS_paths[s,1:T]);
		# This rolling procedure will go all the way until the hurricane gets into the absorbing state of dissipating 
		for t_roll=2:(absorbingT-1)
			# roll up to t = absorbingT-1
			master, x, f, θ, subproblem, y2, xCons, dCons, rCons = RH_2SSP_define_models(t_roll,x_init);
			LB_Roll, UB_Roll, xval_Roll, fval_Roll, θval_Roll = RH_2SSP_solve_roll(s,t_roll,master,subproblem,x,f,θ,y2,xCons,dCons,rCons);
			objGo = UB_Roll;
			let objNoGo = 0.0; # objNoGo: expected obj if wait to implement a two-stage plan in the next stage
				curr_state = allscen[s,t_roll];
				# sample all the next-stage states and see the corresponding expected objval if implementing a two-stage plan then 
				for kk = 1:K
					if P_joint[curr_state, kk] > 1e-8
						# plan with starting state of kk, initial condition of x_init, and t_roll
						master_nextState, x_nextState, f_nextState, θ_nextState, subproblem_nextState, y2_nextState, xCons_nextState, dCons_nextState, rCons_nextState = RH_2SSP_define_models(t_roll+1,x_init);
						LB_nextState, UB_nextState, xval_nextState, fval_nextState, θval_nextState = RH_2SSP_solve_roll2(kk,t_roll+1,master_nextState,subproblem_nextState,x_nextState,f_nextState,θ_nextState,y2_nextState,xCons_nextState,dCons_nextState,rCons_nextState);
						objNoGo += P_joint[curr_state,kk]*UB_nextState;
					end
				end
				if objGo < objNoGo
					Go = true;
				end
			end
			if Go == true
				RH_2SSP_update_RHS(absorbingT,OS_paths[s,absorbingT],subproblem,xCons,dCons,rCons,xval_Roll,fval_Roll,y2,t_roll);
				Qtemp, pi1temp, pi2temp, pi3temp, flagtemp = solve_scen_subproblem(subproblem,xCons,dCons,rCons);
				objs_OOS[s] = LB_Roll-sum(θval_Roll[n]*1.0/nbscen for n = 1:nbscen)+Qtemp;
				break;
			else
				continue;
			end
		end
	else
		# This rolling procedure will stop at t = τ
		for t_roll=2:(τ-1)
			# roll up to t = τ-1
			print("t_roll = ", t_roll);
			println(", τ = ", τ);
			master, x, f, θ, subproblem, y2, xCons, dCons, rCons = RH_2SSP_define_models(t_roll,x_init);
			LB_Roll, UB_Roll, xval_Roll, fval_Roll, θval_Roll = RH_2SSP_solve_roll(s,t_roll,master,subproblem,x,f,θ,y2,xCons,dCons,rCons);
			objGo = UB_Roll;
			let objNoGo = 0.0; # objNoGo: expected obj if wait to implement a two-stage plan in the next stage
				curr_state = allscen[s,t_roll];
				# sample all the next-stage states and see the corresponding expected objval if implementing a two-stage plan then 
				for kk = 1:K
					if P_joint[curr_state, kk] > 1e-8
						# plan with starting state of kk, initial condition of x_init, and t_roll
						master_nextState, x_nextState, f_nextState, θ_nextState, subproblem_nextState, y2_nextState, xCons_nextState, dCons_nextState, rCons_nextState = RH_2SSP_define_models(t_roll+1,x_init);
						LB_nextState, UB_nextState, xval_nextState, fval_nextState, θval_nextState = RH_2SSP_solve_roll2(kk,t_roll+1,master_nextState,subproblem_nextState,x_nextState,f_nextState,θ_nextState,y2_nextState,xCons_nextState,dCons_nextState,rCons_nextState);
						objNoGo += P_joint[curr_state,kk]*UB_nextState;
					end
				end
				if objGo < objNoGo
					Go = true;
				end
			end
			if Go == true
				RH_2SSP_update_RHS(τ,OS_paths[s,τ],subproblem,xCons,dCons,rCons,xval_Roll,fval_Roll,y2,t_roll);
				Qtemp, pi1temp, pi2temp, pi3temp, flagtemp = solve_scen_subproblem(subproblem,xCons,dCons,rCons);
				objs_OOS[s] = LB_Roll-sum(θval_Roll[n]*1.0/nbscen for n = 1:nbscen)+Qtemp;
				break;
			else
				if t_roll < (τ-1)
					continue;
				else
					# Will head into the terminal stage with no preparation!
					for i=1:Ni
						set_normalized_rhs(xCons[i],x_init[i]); # xval_Roll[i,2] gets carried over to period τ 
					end
					k_t = OS_paths[s,t_roll+1];
					for j=1:Nj
						if k_t ∉ absorbing_states           
							set_normalized_rhs(dCons[j], SCEN[k_t][j]);
						else
							set_normalized_rhs(dCons[j], 0);
						end
					end

					set_normalized_rhs(rCons,0); #This is 0 since the cost in the future has not been paid yet -- this is rolling horizon

					# Also need to update the coefficients of y2[i,j] variables in the 2nd stage
					for i=1:Ni
						for j=1:Nj
							set_objective_coefficient(subproblem, y2[i,j], ca[i,j,τ]);
						end
					end
					optimize!(subproblem); 
					objs_OOS[s] += objective_value(subproblem);
				end
			end
		end
	end
	print("obj[", s);
	println("] = ", objs_OOS[s]);
end

elapsed_WS = time() - start;

WS_bar = mean(objs_OOS);
WS_std = std(objs_OOS);
WS_low = WS_bar-1.96*WS_std/sqrt(nbOS);
WS_high = WS_bar+1.96*WS_std/sqrt(nbOS);
println("WS....");
println("μ ± 1.96*σ/√NS = ", WS_bar, " ± ", [WS_low,WS_high]);


fname = "./output/benchmark/rolling2SPresults.csv";
df = CSV.read(fname,DataFrame);

results_WS = Matrix(df);
results_WS[inst,1] = 0;
results_WS[inst,2] = WS_bar;
results_WS[inst,3] = WS_bar-WS_low;
results_WS[inst,4] = 0;
results_WS[inst,5] = elapsed_WS;
results_WS[inst,6] = 0;

updf = DataFrame(results_WS, :auto);
CSV.write(fname,updf);


println("############################################################");
println("############################################################");
