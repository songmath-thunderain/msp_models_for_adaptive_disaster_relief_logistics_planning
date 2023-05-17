s = 1;
t_roll = 1;
x_init = x_0;

osfname = "./data/OOS"*string(k_init)*".csv";
OS_paths = Matrix(CSV.read(osfname,DataFrame)); #read the out-of-sample file
#OS_M = Matrix(CSV.read("./data/inOOS.csv",DataFrame))[:,1] #read the second layer OOS [REVISION: no need anymore]

start = time();

#define the model.
master, x, f, θ, subproblem, y2, xCons, dCons, rCons = RH_2SSP_define_models(t_roll,OS_paths[s,t_roll],x_init);

#solve the model.
LB_1stRoll, UB_1stRoll, xval_1stRoll, fval_1stRoll, θval_1stRoll = RH_2SSP_solve_roll(OS_paths[s,t_roll],t_roll,master,subproblem,x,f,θ,y2,xCons,dCons,rCons);

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
	absorbingT = -1;
	if dissipate_option == 1
		absorbingT = findfirst(x -> (S[x][3] == Nc || S[x][1] == 1), OS_paths[s,1:T]);
	else
		absorbingT = findfirst(x -> (S[x][3] == Nc), OS_paths[s,1:T]);
	end
	x_init = deepcopy(xval_1stRoll[:,1]);
	if dissipate_option == 1 && S[absorbingT][1] == 1
		# This rolling procedure will go all the way until the hurricane gets into the absorbing state of dissipating 
		for t_roll=2:(absorbingT-1)
			# roll up to t = absorbingT-1
			#define the the model.
			master, x, f, θ, subproblem, y2, xCons, dCons, rCons = RH_2SSP_define_models(t_roll,OS_paths[s,t_roll],x_init);
			#solve the model.
			LB_Roll, UB_Roll, xval_Roll, fval_Roll, θval_Roll = RH_2SSP_solve_roll(OS_paths[s,t_roll],t_roll,master,subproblem,x,f,θ,y2,xCons,dCons,rCons);

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
				objs_RH2SSP[s,t_roll] = 0;
				# Just salvage all the x_init

				# Regardless of what the absorbing_option is, we won't do anything but to salvage everything
				for i=1:Ni
					objs_RH2SSP[s,t_roll] += xval_Roll[i,1]*q;
				end
			end
		end
	else
		# This rolling procedure will stop at t = absorbingT
		for t_roll=2:(absorbingT-1)
			# roll up to t = absorbingT-1
			#define the the model.
			master, x, f, θ, subproblem, y2, xCons, dCons, rCons = RH_2SSP_define_models(t_roll,OS_paths[s,t_roll],x_init);
			#solve the model.
			LB_Roll, UB_Roll, xval_Roll, fval_Roll, θval_Roll = RH_2SSP_solve_roll(OS_paths[s,t_roll],t_roll,master,subproblem,x,f,θ,y2,xCons,dCons,rCons);
	
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
			
			if t_roll == (absorbingT-1)
				# Now we get the realization, do the recourse now and finish the rolling procedure
				t_roll = t_roll + 1;
				objs_RH2SSP[s,t_roll] = 0;

				# Approach #2: based on xvals_Roll[:,1], optimize all the operations together with full information
				if absorbing_option == 0
					# do not allow additional MDP/SP operation
					for i=1:Ni
						set_normalized_rhs(xCons[i],xval_Roll[i,1]); 
					end
					k_t = OS_paths[s,t_roll];
					for j=1:Nj
						if S[k_t][1] != 1           
							set_normalized_rhs(dCons[j], SCEN[k_t][j]);
						else
							set_normalized_rhs(dCons[j], 0);
						end
					end

					set_normalized_rhs(rCons,0); #This is 0 since the cost in the future has not been paid yet -- this is rolling horizon

					# Also need to update the coefficients of y2[i,j] variables in the 2nd stage
					for i=1:Ni
						for j=1:Nj
							set_objective_coefficient(subproblem, y2[i,j], ca[i,j,absorbingT]);
						end
					end
					optimize!(subproblem); 
       				if termination_status(subproblem) != MOI.OPTIMAL
						println("status = ", termination_status(subproblem));
						exit(0);
					else
						objs_RH2SSP[s,t_roll] += objective_value(subproblem);
					end	
				else
					# Solve a terminal stage problem, just as the FA/MSP version
					m_term, x_term, f_term, y_term, z_term, v_term, dCons_term = terminal_model(t_roll,xval_Roll[:,1]);
					k_t = OS_paths[s,t_roll];
					for j=1:Nj
						if S[k_t][1] != 1           
							set_normalized_rhs(dCons_term[j], SCEN[k_t][j]);
						else
							set_normalized_rhs(dCons_term[j], 0);
						end
					end
					optimize!(m_term); 
       				if termination_status(m_term) != MOI.OPTIMAL
						println("status = ", termination_status(m_term));
						exit(0);
					else
						objs_RH2SSP[s,t_roll] += objective_value(m_term);
					end
				end

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
