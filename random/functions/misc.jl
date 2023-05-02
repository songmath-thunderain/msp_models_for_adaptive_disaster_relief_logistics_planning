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
	tIter = 1;
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
		#println("tempList = ", tempList);
		push!(nodeLists,tempList);			
	end
	return nodeLists;
end
