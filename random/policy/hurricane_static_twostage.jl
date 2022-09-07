#max_iter = 100;
s = 1;
t_roll = 1;
x_init = x_0;
#define the the model.
master, x, f, θ, subproblem, y2, xCons, dCons, rCons = RH_2SSP_define_models(t_roll,x_init);

#solve the model.
start=time();
LB_st2SSP, UB_st2SSP, xval_st2SSP, fval_st2SSP, θval_st2SSP = RH_2SSP_solve_roll(s,t_roll,master,subproblem,x,f,θ,y2,xCons,dCons,rCons);
elapsed = time() - start;

f1cost = LB_st2SSP-θval_st2SSP;

#eval the model.
start=time();

OS_paths = Matrix(CSV.read("./data/OOS.csv",DataFrame)); #read the out-of-sample file
#OS_M = Matrix(CSV.read("./data/inOOS.csv",DataFrame))[:,1] #read the second layer OOS [REVISION: no need any more]
    
objs_st2SSP = fill(f1cost,nbOS);
Q = zeros(nbOS); #list for all the optimal values
pi1 = Array{Any,1}(undef,nbOS); 
pi2 = Array{Any,1}(undef,nbOS); 
pi3 = zeros(nbOS);

for s=1:nbOS
	#identify the period where the hurricane makes landfall 
	τ = findfirst(x -> S[x][3] == Nc-1 && x ∉ absorbing_states, OS_paths[s,1:T]);
	
	#update the RHS
	if τ === nothing
		#RH_2SSP_update_RHS(τ,OS_paths[s,end],OS_M[s],nbstages1,ConsFB,dCons,rCons,xval,fval,t_roll) [REVISION]
		RH_2SSP_update_RHS(τ,OS_paths[s,end],subproblem,xCons,dCons,rCons,xval_st2SSP,fval_st2SSP,y2,t_roll)
	else
		#RH_2SSP_update_RHS(τ,OS_paths[s,τ],OS_M[s],nbstages1,ConsFB,dCons,rCons,xval,fval,t_roll) [REVISION]
		RH_2SSP_update_RHS(τ,OS_paths[s,τ],subproblem,xCons,dCons,rCons,xval_st2SSP,fval_st2SSP,y2,t_roll)
	end
	
	#solve the subproblem and store the dual information
	flag = solve_scen_subproblem(Q,pi1,pi2,pi3,s,subproblem,xCons,dCons,rCons)
	objs_st2SSP[s] = objs_st2SSP[s] + Q[s];
end
  
st2SSP_bar = mean(objs_st2SSP);
st2SSP_std = std(objs_st2SSP);
st2SSP_low = st2SSP_bar-1.96*st2SSP_std/sqrt(nbOS);
st2SSP_high = st2SSP_bar+1.96*st2SSP_std/sqrt(nbOS);
println("μ ± 1.96*σ/√NS = ", st2SSP_bar, " ± ", [st2SSP_low,st2SSP_high]);

elapsed2 = time() - start;

fname = "./output/benchmark/static2SPresults.csv"
df = CSV.read(fname,DataFrame);

results_st2SSP = Matrix(df);
results_st2SSP[inst,1] = LB_st2SSP
results_st2SSP[inst,2] = st2SSP_bar
results_st2SSP[inst,3] = st2SSP_bar-st2SSP_low
results_st2SSP[inst,4] = elapsed
results_st2SSP[inst,5] = elapsed2
results_st2SSP[inst,6] = 0

updf = DataFrame(results_st2SSP, :auto);
CSV.write(fname,updf)
