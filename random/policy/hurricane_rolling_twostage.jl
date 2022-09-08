s = 1;
t_roll = 1;
x_init = x_0;

start = time();

#define the model.
master, x, f, θ, subproblem, y2, xCons, dCons, rCons = RH_2SSP_define_models(t_roll,x_init);

#solve the model.
LB_1stRoll, UB_1stRoll, xval_1stRoll, fval_1stRoll, θval_1stRoll = RH_2SSP_solve_roll(s,t_roll,master,subproblem,x,f,θ,y2,xCons,dCons,rCons);

OS_paths = Matrix(CSV.read("./data/OOS.csv",DataFrame)); #read the out-of-sample file
#OS_M = Matrix(CSV.read("./data/inOOS.csv",DataFrame))[:,1] #read the second layer OOS [REVISION: no need anymore]

f1cost = LB_1stRoll-θval_1stRoll;

objs_RH2SSP = zeros(nbOS,T);
objs_RH2SSP[:,1] .= f1cost;

for s=1:nbOS
	τ = findfirst(x -> S[x][3] == Nc-1 && x ∉ absorbing_states, OS_paths[s,1:T]);
	x_init = deepcopy(xval_1stRoll[:,1]);
	# Here we may seet aside the special case when the hurricane dissipates [TBD]
	if τ !== nothing
		# This rolling procedure will stop at t = τ
		for t_roll=2:(τ-1)
			# roll up to t = τ-1
			#define the the model.
			master, x, f, θ, subproblem, y2, xCons, dCons, rCons = RH_2SSP_define_models(t_roll,x_init);
			
			#solve the model.
			LB_Roll, UB_Roll, xval_Roll, fval_Roll, θval_Roll = RH_2SSP_solve_roll(s,t_roll,master,subproblem,x,f,θ,y2,xCons,dCons,rCons);

			#implement xₜ, pay cost, and pass xₜ to new t+1.
			x_init = deepcopy(xval_Roll[:,1]);
			objs_RH2SSP[s,t_roll] = LB_Roll - θval_Roll;
			if t_roll == (τ-1)
				# Now we get the realization, do the recourse now and finish the rolling procedure
				t_roll = t_roll + 1;
				nbstages1 = T-t_roll+1;
				for i=1:Ni
					set_normalized_rhs(xCons[i],xval_Roll[i,2]); # note that here index 2 is a relative index, the reason why it is 2 is because xval_Roll comes from the previous roll
				end
				k_t = OS_paths[s,t_roll];
				for j=1:Nj
					if k_t ∉ absorbing_states           
						set_normalized_rhs(dCons[j], SCEN[k_t][j]);
					else
						set_normalized_rhs(dCons[j], 0);
					end
				end

				set_normalized_rhs(rCons,
								-sum(sum(sum(cb[i,ii,t_roll+t-1]*fval_Roll[i,ii,t] for ii=1:Ni) for i=1:N0)
								+sum(ch[i,t_roll+t-1]*xval_Roll[i,t] for i=1:Ni)  
								+sum(fval_Roll[N0,i,t] for i=1:Ni)*h[t_roll+t-1] for t = (τ+2-t_roll):nbstages1) 
								);
				# Also need to update the coefficients of y2[i,j] variables in the 2nd stage
				for i=1:Ni
					for j=1:Nj
						set_objective_coefficient(subproblem, y2[i,j], ca[i,j,τ]);
					end
				end
				optimize!(subproblem); 
      			objs_RH2SSP[s,t_roll] = objective_value(subproblem);
			end
		end
	else
		tt = findfirst(x -> S[x][1] == 1 && x ∈ absorbing_states, OS_paths[s,1:T]);
		# This rolling procedure will go all the way to the end
		for t_roll=2:(T-1)
			# roll up to t = T-1
			#define the the model.
			master, x, f, θ, subproblem, y2, xCons, dCons, rCons = RH_2SSP_define_models(t_roll,x_init);
			#solve the model.
			LB_Roll, UB_Roll, xval_Roll, fval_Roll, θval_Roll = RH_2SSP_solve_roll(s,t_roll,master,subproblem,x,f,θ,y2,xCons,dCons,rCons);
			#implement xₜ, pay cost, and pass xₜ to new t+1.
			x_init = deepcopy(xval_Roll[:,1]);
			objs_RH2SSP[s,t_roll] = LB_Roll - θval_Roll;
			if t_roll == (T-1)
				# Now we get the realization, do the recourse now and finish the rolling procedure
				t_roll = t_roll + 1;
				objs_RH2SSP[s,t_roll] = θval_Roll; # all scenarios should give the same result
			end
		end
	end
end

elapsed_RH2SSP = time() - start;

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
