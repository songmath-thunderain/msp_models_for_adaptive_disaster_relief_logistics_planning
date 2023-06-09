#Fully adaptive
#Defne the model 
m_fa, x_fa, f_fa, y_fa, z_fa, v_fa, ϴ_fa, dCons_fa, FB1Cons_fa, FB2Cons_fa = define_models();


#train the model
LB, train_time, iter = train_models_offline();

println("***********************************************")
println("***********************************************")
println("finished training!")
println("LB = ", LB[end])
println("training time = ", train_time)
println("number of iterations = ", iter)

#evaluate the model 
start=time();
osfname = "./data/OOS"*string(k_init)*".csv";
OS_paths = Matrix(CSV.read(osfname,DataFrame)); #read the out-of-sample file
#OS_M = Matrix(CSV.read("./data/inOOS.csv",DataFrame))[:,1] #read the second layer OOS
objs_fa = zeros(nbOS,T);
procurmnt_all = zeros(nbOS,T);
xval_fa = Array{Any,2}(undef,nbOS,T); fval_fa = Array{Any,2}(undef,nbOS,T);
yval_fa = Array{Any,2}(undef,nbOS,T); zval_fa = Array{Any,2}(undef,nbOS,T); vval_fa = Array{Any,2}(undef,nbOS,T);
procurmnt_amount = zeros(T); 
procurmnt_percentage = zeros(T); 
procurmnt_posExpect = zeros(T); 
flow_amount = zeros(T);
invAmount = zeros(nbOS);
salvageAmount = zeros(nbOS);
penaltyAmount = zeros(nbOS);
for s=1:nbOS
    xval = zeros(Ni,T);
    for t=1:T
        #the state is known in the first stage; if not sample a new state k 
        k_t = OS_paths[s,t];
		# we do not have this second layer now [REVISION]
		#m = OS_M[s]; # realization from OS path corresponding to layer 2
    
        if t > 1
            MSP_fa_update_RHS(k_t,t,xval);
        end
        #solve the model
        optimize!(m_fa[t,k_t])

        #check the status 
        status = termination_status(m_fa[t,k_t]);
        if status != MOI.OPTIMAL
            println(" in evaluation")
            println("Model in stage =", t, " and state = ", k_t, ", in forward pass is ", status)
            return 
        else
            #collect values
            xval_fa[s,t] = value.(x_fa[t,k_t]); 
            xval[:,t] = xval_fa[s,t];
            fval_fa[s,t] = value.(f_fa[t,k_t]);                 
            yval_fa[s,t] = value.(y_fa[t,k_t]);
            zval_fa[s,t] = value.(z_fa[t,k_t]);
            vval_fa[s,t] = value.(v_fa[t,k_t]);
            objs_fa[s,t] = objective_value(m_fa[t,k_t])- value(ϴ_fa[t,k_t]);
            if (k_t in absorbing_states) == false
                procurmnt_all[s,t] = sum(fval_fa[s,t][N0,i] for i=1:Ni);
                procurmnt_amount[t] += (sum(fval_fa[s,t][N0,i] for i=1:Ni))/nbOS;
                flow_amount[t] += sum(sum(fval_fa[s,t][i,ii] for i = 1:Ni) for ii = 1:Ni)/nbOS;
            else
                invAmount[s] = sum(xval_fa[s,t-1][i] for i=1:Ni);
            end
            salvageAmount[s] += sum(vval_fa[s,t][i] for i=1:Ni);
            penaltyAmount[s] += sum(zval_fa[s,t][j] for j=1:Nj);
        end
        if k_t in absorbing_states
            break;
        end
    end        
end

for t=1:T
	count = 0;
	totalPos = 0;
	for s=1:nbOS
		if procurmnt_all[s,t] > 1e-2
			count += 1;
			totalPos += procurmnt_all[s,t];
		end
	end
	procurmnt_percentage[t] = count*1.0/nbOS;
	if count > 0
		procurmnt_posExpect[t] = totalPos*1.0/count;
	end
end

println("procurement amount = ", procurmnt_amount);
println("flow amount = ", flow_amount);

avgInvAmount = sum(invAmount)*1.0/nbOS;
avgSalvageAmount = sum(salvageAmount)*1.0/nbOS;
avgPenaltyAmount = sum(penaltyAmount)*1.0/nbOS;


println("avgInvAmount = ", avgInvAmount);
println("avgSalvageAmount = ", avgSalvageAmount);
println("avgPenaltyAmount = ", avgPenaltyAmount);

fname = "./output/FA-sensitivity.csv"
df = CSV.read(fname,DataFrame);
results_fa = Matrix(df);
results_fa[inst,1:T] = procurmnt_amount
results_fa[inst,(T+1):(2*T)] = procurmnt_percentage
results_fa[inst,(2*T+1):(3*T)] = procurmnt_posExpect
results_fa[inst,(3*T+1):(4*T)] = flow_amount
results_fa[inst,4*T+2] = avgInvAmount
results_fa[inst,4*T+3] = avgSalvageAmount
results_fa[inst,4*T+4] = avgPenaltyAmount
results_fa[inst,end] = inst; 
updf = DataFrame(results_fa, :auto);
CSV.write(fname,updf)
