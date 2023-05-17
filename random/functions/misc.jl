#miscellaneous functions 

#sample a markovian state
function MC_sample(current_state)
    states = collect(1:1:K);
    k = sample(states, Weights(P_joint[current_state,:]));
    return k
end


#training termination check
function termination_check(iter,relative_gap,LB,start,cutviol_iter)
    flag = 0;
    Elapsed = time() - start;
    if iter > max_iter
        flag = 1
        println("max iteration is reached")
    elseif Elapsed > time_limit
        flag = 2
        println("time limit is reached")
    elseif cutviol_iter > cutviol_maxiter
        flag = 3
        println("cut violation is reached")
    else
        if iter > stall
            relative_gap = (LB[iter]-LB[iter-stall])/max(1e-10,abs(LB[iter-stall]));
            if relative_gap < ϵ
                flag = 4
                println("the LB is not making significant progress");
            end
        end
    end
    return flag, Elapsed
end


# function to create an out-of-sample for given initial state k_init
function create_OSpaths(k_init)
	osfname = "./data/OOS"*string(k_init)*".csv";
    OS_paths = Matrix(CSV.read(osfname,DataFrame)); #read the out-of-sample file
    rr = size(OS_paths)[1];
    cc = size(OS_paths)[2];
    OS_paths = fill(k_init,rr,cc);
    for s=1:rr
        for t=2:cc
            OS_paths[s,t] = MC_sample(OS_paths[s,t-1]) 
        end
    end
    df = DataFrame(OS_paths, :auto);
    CSV.write(osfname,df)
end

#function to save lp files
function save_lp(lp,name)
    open("./output/lp/$name.lp", "w") do f
        print(f, lp)
    end
end

#function save an n×m matrix as an csv file
function save_csv(m_data, m_colname, f_dir, m_fname)
    fname = f_dir*m_fname*".csv"
    df = DataFrame(m_data, m_colname);
    CSV.write(fname,df)
end

# function that create a complete set of reachable nodes over time, starting from the initial state k_init
function createNodes(k_init)
	nodeLists = Array{Int}[];
	push!(nodeLists,[k_init]);
	stopFlag = false;
	while stopFlag == false
		tempList = [];
		stopFlag = true; # if there is at least one state that is absorbing, turn the stopFlag back to false
		for k=1:K
			for kk in last(nodeLists)
				#@printf("kk = %d, k = %d, P_joint = %f, smallestTransProb = %f", kk, k, P_joint[kk,k],smallestTransProb);
				if (kk in absorbing_states) == false && P_joint[kk,k] > smallestTransProb
					push!(tempList,k);
					if (k in absorbing_states) == false
						stopFlag = false;
					end
					break;
				end
			end
		end
		push!(nodeLists,tempList);			
	end
	return nodeLists;
end

# function that creates a list of scenarios, along with the probability of occurrence, for each transient state node in the nodeList (set of reachable nodes from the initial state k_init)
function createNodeScens(k_init, nodeLists)
	# Note: there might be repetitions in the nodeLists, e.g., the same Hurricane state can be reached at different times
	nodeScenList = Dict();
	nodeScenWeights = Dict();

	for t = (T-1):-1:1
		for k = 1:length(nodeLists[t])
			if (nodeLists[t][k] in absorbing_states) == false
				nodeScenList[t,nodeLists[t][k]] = [];
				nodeScenWeights[t,nodeLists[t][k]] = [];
				for kk = 1:length(nodeLists[t+1])
					if P_joint[nodeLists[t][k],nodeLists[t+1][kk]] > smallestTransProb
						if nodeLists[t+1][kk] in absorbing_states
							# absorbing states, directly append 
							push!(nodeScenList[t,nodeLists[t][k]],[t+1,nodeLists[t+1][kk]]);
							push!(nodeScenWeights[t,nodeLists[t][k]], P_joint[nodeLists[t][k],nodeLists[t+1][kk]]);
						else
							# transient states, append the corresponding scenlist and weights
							for j = 1:length(nodeScenList[t+1,nodeLists[t+1][kk]])
								if (nodeScenList[t+1,nodeLists[t+1][kk]][j] in nodeScenList[t,nodeLists[t][k]]) == false
									# Not in the scenario list, so go ahead and add it
									push!(nodeScenList[t,nodeLists[t][k]],nodeScenList[t+1,nodeLists[t+1][kk]][j]);
									push!(nodeScenWeights[t,nodeLists[t][k]],P_joint[nodeLists[t][k],nodeLists[t+1][kk]]*nodeScenWeights[t+1,nodeLists[t+1][kk]][j]);
								else
									# in the scenario list, increment the probability
									ind = findfirst(x->x==nodeScenList[t+1,nodeLists[t+1][kk]][j],nodeScenList[t,nodeLists[t][k]]);
									nodeScenWeights[t,nodeLists[t][k]][ind] += P_joint[nodeLists[t][k],nodeLists[t+1][kk]]*nodeScenWeights[t+1,nodeLists[t+1][kk]][j];
								end
							end
						end
					end
				end	
				if abs(sum(nodeScenWeights[t,nodeLists[t][k]])-1) > 1e-6
					println("Wrong!");
					exit(0);
				end
				#println("node[", nodeLists[t][k], "]'s scenario # = ", length(nodeScenList[nodeLists[t][k]]));
				#if t == 1
				#	println("scenlist = ", nodeScenList[nodeLists[t][k]]);
				#	println("prob = ", nodeScenWeights[nodeLists[t][k]]);
				#end
			end
		end
	end

	return nodeScenList, nodeScenWeights;
end
