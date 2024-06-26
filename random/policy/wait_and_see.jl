start = time();

###############################################################
###############################################################
# Offline policy construction: for each node in the nodeLists, going backwards in time, choose Go/No-go and record the best (better) expected objective value
# The expected objective value of No-go nodes: a simple calculation for an absorbing state, and an expectation over all reachable nodes for a transient state

# initialize -- just a placeholder
decisionNodes = deepcopy(nodeLists);
objvalNodes = Array{Float64}[];
for t = 1:T
	push!(objvalNodes,zeros(length(nodeLists[t])));
end	
solutionNodes = Dict();

for t_roll = 1:T
	for k = 1:length(nodeLists[t_roll])
		if (nodeLists[t_roll][k] in absorbing_states) == false
			# transient state: solving a 2SSP
			x_init = deepcopy(x_0);
			decisionNodes[t_roll][k] = 1; # This is only temporary, need to do a round of backward cleanup
			#define the model.
			master, x, f, θ, subproblem, y2, xCons, dCons, rCons = RH_2SSP_define_models(t_roll,nodeLists[t_roll][k],x_init);
			#solve the model.
			LB_Roll, UB_Roll, xval_Roll, fval_Roll, θval_Roll = RH_2SSP_solve_roll(nodeLists[t_roll][k],t_roll,master,subproblem,x,f,θ,y2,xCons,dCons,rCons);
			objvalNodes[t_roll][k] = UB_Roll; # objGo: expected objval if implementing a two-stage plan now. This is also only temporary, need to do a round of backward cleanup
			solutionNodes[t_roll,k] = [xval_Roll, fval_Roll];
		else
			# absorbing state, we can determine pretty easily if it is Go vs. No-Go
			if S[nodeLists[t_roll][k]][1] == 1 && dissipate_option == 1
				# Hurricane dissipates
				objvalNodes[t_roll][k] = 0;
				decisionNodes[t_roll][k] = 0;
			else
				# Hurricane makes landfall with intensity, need to decide to do nothing or do some last-minute operation (although perhaps costly)
				costNoGo = p*sum(SCEN[nodeLists[t_roll][k]][j] for j = 1:Nj);
				if absorbing_option == 0
					objvalNodes[t_roll][k] = costNoGo;
					decisionNodes[t_roll][k] = 0;
				else
					#define terminal stage optimality model

					m_term, x_term, f_term, y_term, z_term, v_term, dCons_term = terminal_model(t_roll,x_0);

					for j=1:Nj
						set_normalized_rhs(dCons_term[j], SCEN[nodeLists[t_roll][k]][j]);
					end

					optimize!(m_term) #solve the model
					status_subproblem = termination_status(m_term); #check the status 
					costGo = 0;
					if status_subproblem != MOI.OPTIMAL
						println("status_subproblem = ", status_subproblem);
						exit(0);
					else
						costGo = objective_value(m_term);
					end
					if costNoGo < costGo
						objvalNodes[t_roll][k] = costNoGo;
						decisionNodes[t_roll][k] = 0;
					else
						objvalNodes[t_roll][k] = costGo;
						decisionNodes[t_roll][k] = 1;
					end
				end
			end
		end
	end
end


###############################################################
###############################################################
# Now let's do a round of backward correction!
for t=(T-1):-1:1
	# Starting from T-1
	for k = 1:length(nodeLists[t])
		if (nodeLists[t][k] in absorbing_states) == false
			# we've computed the Go version for these nodes above, now let's see the no-go version
			costNoGo = 0;
			for kk = 1:length(nodeLists[t+1])
				if P_joint[nodeLists[t][k],nodeLists[t+1][kk]] > smallestTransProb
					costNoGo = costNoGo + P_joint[nodeLists[t][k],nodeLists[t+1][kk]]*objvalNodes[t+1][kk];
				end
			end
			if costNoGo < objvalNodes[t][k]
				decisionNodes[t][k] = 0;
				objvalNodes[t][k] = costNoGo;
			end
		end
	end
end


#for t = 1:T
#	println("objvalNodes[", t, "] = ", objvalNodes[t]);
#end

###############################################################
###############################################################
# Start evaluating policies on the decision tree

println("Construction is done....Now we do evaluation...");

osfname = "./data/OOS"*string(k_init)*".csv";
OS_paths = Matrix(CSV.read(osfname,DataFrame)); #read the out-of-sample file
objs_OOS = zeros(nbOS);

for s=1:nbOS
	print("OOS[", s, "] = (");
	for t = 1:T
		ind = findfirst(x->x==OS_paths[s,t], nodeLists[t]);
		print("t = ", t, ", ind = ", ind);
		if OS_paths[s,t] in absorbing_states
			# if absorbing, just take whatever that is the best, which has been computed above
			objs_OOS[s] = objvalNodes[t][ind];
			print("absorbed! obj = ", objs_OOS[s], "\n");
			break;
		else
			if decisionNodes[t][ind] == 0
				# Decision is No-go!
				#print("No-go! ");
				continue;
			else
				# Decision is Go!
				absorbingT = -1;
				if dissipate_option == 1
					absorbingT = findfirst(x -> (S[x][3] == Nc || S[x][1] == 1), OS_paths[s,:]);
				else
					absorbingT = findfirst(x -> (S[x][3] == Nc), OS_paths[s,:]);
				end
				#define second stage (subproblem) optimality model
   				subproblem, y, xCons, dCons, rCons = RH_2SSP_second_stage();
				if absorbing_option == 0
					for i=1:Ni
						set_normalized_rhs(xCons[i],solutionNodes[t,ind][1][i,absorbingT-t]);
					end
				else	
					for i=1:Ni
						set_normalized_rhs(xCons[i],solutionNodes[t,ind][1][i,absorbingT-t+1]);
					end
				end
				for j=1:Nj
					set_normalized_rhs(dCons[j], SCEN[OS_paths[s,absorbingT]][j]);
				end

				for i=1:Ni
					for j=1:Nj
						set_objective_coefficient(subproblem, y[i,j], ca[i,j,absorbingT]);
					end
				end

				optimize!(subproblem) #solve the model
   				status_subproblem = termination_status(subproblem); #check the status 
   		 		if status_subproblem != MOI.OPTIMAL
					println("status_subproblem = ", status_subproblem);
					exit(0);
    			else
        			objs_OOS[s] = objective_value(subproblem);
					print("first obj = ", objs_OOS[s]);
					if absorbing_option == 0
						for tt = 1:(absorbingT-t)
							objs_OOS[s] = objs_OOS[s] + (sum(sum(cb[i,ii,t+tt-1]*solutionNodes[t,ind][2][i,ii,tt] for ii=1:Ni) for i=1:N0)+sum(ch[i,t+tt-1]*solutionNodes[t,ind][1][i,tt] for i=1:Ni)+sum(solutionNodes[t,ind][2][N0,i,tt] for i=1:Ni)*h[t+tt-1]);
						end
					else
						for tt = 1:(absorbingT+1-t)
							objs_OOS[s] = objs_OOS[s] + (sum(sum(cb[i,ii,t+tt-1]*solutionNodes[t,ind][2][i,ii,tt] for ii=1:Ni) for i=1:N0)+sum(ch[i,t+tt-1]*solutionNodes[t,ind][1][i,tt] for i=1:Ni)+sum(solutionNodes[t,ind][2][N0,i,tt] for i=1:Ni)*h[t+tt-1]);
						end
					end
				end
				print("Go! obj = ", objs_OOS[s], "\n");
				break;
			end
		end
	end
end


elapsed_WS = time() - start;

WS_bar = mean(objs_OOS);
WS_std = std(objs_OOS);
WS_low = WS_bar-1.96*WS_std/sqrt(nbOS);
WS_high = WS_bar+1.96*WS_std/sqrt(nbOS);
println("WS....");
println("μ ± 1.96*σ/√NS = ", WS_bar, " ± ", [WS_low,WS_high]);

fname = "./output/benchmark/wait-and-see.csv";
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
