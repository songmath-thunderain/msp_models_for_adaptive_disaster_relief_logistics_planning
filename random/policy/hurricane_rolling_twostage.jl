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

objs_RH2SSP = zeros(nbOS,T);

for s=1:nbOS
	for i = 1:N0
		for ii = 1:Ni
			objs_RH2SSP[s,1] += cb[i,ii,1]*fval_1stRoll[i,ii,1];
		end
	end
	for i = 1:Ni
		objs_RH2SSP[s,1] += (ch[i,1]*xval_1stRoll[i,1] + h[1]*fval_1stRoll[N0,i,1]);
	end
end

for s=1:nbOS
	τ = findfirst(x -> S[x][3] == Nc-1 && x ∉ absorbing_states, OS_paths[s,1:T]);
	x_init = deepcopy(xval_1stRoll[:,1]);
	if τ === nothing
		absorbingT = findfirst(x -> S[x][1] == 1, OS_paths[s,1:T]);
		# This rolling procedure will go all the way until the hurricane gets into the absorbing state of dissipating 
		for t_roll=2:(absorbingT-1)
			# roll up to t = absorbingT-1
			#define the the model.
			master, x, f, θ, subproblem, y2, xCons, dCons, rCons = RH_2SSP_define_models(t_roll,x_init);
			#solve the model.
			LB_Roll, UB_Roll, xval_Roll, fval_Roll, θval_Roll = RH_2SSP_solve_roll(s,t_roll,master,subproblem,x,f,θ,y2,xCons,dCons,rCons);
#=
			println("roll #", t_roll);
			for t=1:(T-t_roll+1)
				print("xval[", t);
				println("] = ", xval_Roll[:,t]);
			end
=#
			#implement xₜ, pay cost, and pass xₜ to new t+1.
			x_init = deepcopy(xval_Roll[:,1]);
			objs_RH2SSP[s,t_roll] = 0;
			for i = 1:N0
				for ii = 1:Ni
					objs_RH2SSP[s,t_roll] += cb[i,ii,t_roll]*fval_Roll[i,ii,1];
				end
			end
			for i = 1:Ni
				objs_RH2SSP[s,t_roll] += (ch[i,t_roll]*xval_Roll[i,1] + h[t_roll]*fval_Roll[N0,i,1]);
			end
			if t_roll == (absorbingT-1)
				# Now we get the realization, do the recourse now and finish the rolling procedure
				t_roll = t_roll + 1;
				# Just salvage all the x_init
				objs_RH2SSP[s,t_roll] = 0;
				for i = 1:N0
					for ii = 1:Ni
						objs_RH2SSP[s,t_roll] += cb[i,ii,t_roll]*fval_Roll[i,ii,2];
					end
				end
				for i = 1:Ni
					objs_RH2SSP[s,t_roll] += (ch[i,t_roll]*xval_Roll[i,2] + h[t_roll]*fval_Roll[N0,i,2]);
				end
				for i=1:Ni
					objs_RH2SSP[s,t_roll] += xval_Roll[i,2]*q;
				end
			end
		end
	else
		# This rolling procedure will stop at t = τ
		for t_roll=2:(τ-1)
			# roll up to t = τ-1
			#define the the model.
			master, x, f, θ, subproblem, y2, xCons, dCons, rCons = RH_2SSP_define_models(t_roll,x_init);
			#solve the model.
			LB_Roll, UB_Roll, xval_Roll, fval_Roll, θval_Roll = RH_2SSP_solve_roll(s,t_roll,master,subproblem,x,f,θ,y2,xCons,dCons,rCons);
#=
			println("roll #", t_roll);
			for t=1:(T-t_roll+1)
				print("xval[", t);
				println("] = ", xval_Roll[:,t]);
			end
=#
			#implement xₜ, pay cost, and pass xₜ to new t+1.
			x_init = deepcopy(xval_Roll[:,1]);
			# note that we should only store the cost incurred at the current period: the first-stage cost includes other periods!
			objs_RH2SSP[s,t_roll] = 0;
			for i = 1:N0
				for ii = 1:Ni
					objs_RH2SSP[s,t_roll] += cb[i,ii,t_roll]*fval_Roll[i,ii,1];
				end
			end
			for i = 1:Ni
				objs_RH2SSP[s,t_roll] += (ch[i,t_roll]*xval_Roll[i,1] + h[t_roll]*fval_Roll[N0,i,1]);
			end
			
			if t_roll == (τ-1)
				# Now we get the realization, do the recourse now and finish the rolling procedure
				t_roll = t_roll + 1;
				objs_RH2SSP[s,t_roll] = 0;
				for i = 1:N0
					for ii = 1:Ni
						objs_RH2SSP[s,t_roll] += cb[i,ii,t_roll]*fval_Roll[i,ii,2];
					end
				end
				for i = 1:Ni
					objs_RH2SSP[s,t_roll] += (ch[i,t_roll]*xval_Roll[i,2] + h[t_roll]*fval_Roll[N0,i,2]);
				end
				for i=1:Ni
					set_normalized_rhs(xCons[i],xval_Roll[i,2]); # xval_Roll[i,2] gets carried over to period τ 
				end
				k_t = OS_paths[s,t_roll];
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
      			objs_RH2SSP[s,t_roll] += objective_value(subproblem);
				#println("last roll value = ", objs_RH2SSP[s,t_roll]);
			end
		end
	end
	print("obj[", s);
	println("] = ", sum(objs_RH2SSP[s,t] for t=1:T));
end

elapsed_RH2SSP = time() - start;

RH2SSP_bar = mean(sum(objs_RH2SSP[:,t] for t=1:T));
RH2SSP_std = std(sum(objs_RH2SSP[:,t] for t=1:T));
RH2SSP_low = RH2SSP_bar-1.96*RH2SSP_std/sqrt(nbOS);
RH2SSP_high = RH2SSP_bar+1.96*RH2SSP_std/sqrt(nbOS);
println("RH2SSP....");
println("μ ± 1.96*σ/√NS = ", RH2SSP_bar, " ± ", [RH2SSP_low,RH2SSP_high]);


fname = "./output/benchmark/rolling2SPresults.csv";
df = CSV.read(fname,DataFrame);

results_RH2SSP = Matrix(df);
results_RH2SSP[inst,1] = 0;
results_RH2SSP[inst,2] = RH2SSP_bar;
results_RH2SSP[inst,3] = RH2SSP_bar-RH2SSP_low;
results_RH2SSP[inst,4] = 0;
results_RH2SSP[inst,5] = elapsed_RH2SSP;
results_RH2SSP[inst,6] = 0;

updf = DataFrame(results_RH2SSP, :auto);
CSV.write(fname,updf);


println("############################################################");
println("############################################################");
