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

f1cost = LB_st2SSP-sum(θval_st2SSP[n] for n = 1:nbscen)*1.0/nbscen;

println("training LB = ", LB_st2SSP);
println("training UB = ", UB_st2SSP);

#eval the model.
start=time();

osfname = "./data/OOS"*string(k_init)*".csv";
OS_paths = Matrix(CSV.read(osfname,DataFrame)); #read the out-of-sample file

#OS_M = Matrix(CSV.read("./data/inOOS.csv",DataFrame))[:,1] #read the second layer OOS [REVISION: no need any more]
    
objs_st2SSP = fill(f1cost,nbOS);
Q = zeros(nbOS); #list for all the optimal values
pi1 = Array{Any,1}(undef,nbOS); 
pi2 = Array{Any,1}(undef,nbOS); 
pi3 = zeros(nbOS);

for s=1:nbOS
	#identify the period where the hurricane makes landfall 
	absorbingT = findfirst(x -> (S[x][3] == Nc || S[x][1] == 1), OS_paths[s,:]);
	RH_2SSP_update_RHS(absorbingT,OS_paths[s,absorbingT],subproblem,xCons,dCons,rCons,xval_st2SSP,fval_st2SSP,y2,t_roll);
	
	#solve the subproblem and store the dual information
	Q[s], pi1[s], pi2[s], pi3[s], flag = solve_scen_subproblem(subproblem,xCons,dCons,rCons);
	objs_st2SSP[s] = objs_st2SSP[s] + Q[s];
	print("obj[", s);
	println("] = ", objs_st2SSP[s]);
end
  
st2SSP_bar = mean(objs_st2SSP);
st2SSP_std = std(objs_st2SSP);
st2SSP_low = st2SSP_bar-1.96*st2SSP_std/sqrt(nbOS);
st2SSP_high = st2SSP_bar+1.96*st2SSP_std/sqrt(nbOS);
println("static 2SSP....");
println("μ ± 1.96*σ/√NS = ", st2SSP_bar, " ± ", [st2SSP_low,st2SSP_high]);

elapsed2 = time() - start;

fname = "./output/benchmark/static2SPresults.csv"
df = CSV.read(fname,DataFrame);

results_st2SSP = Matrix(df);
results_st2SSP[inst,1] = LB_st2SSP;
results_st2SSP[inst,2] = st2SSP_bar;
results_st2SSP[inst,3] = st2SSP_bar-st2SSP_low;
results_st2SSP[inst,4] = elapsed;
results_st2SSP[inst,5] = elapsed2;
results_st2SSP[inst,6] = 0;

updf = DataFrame(results_st2SSP, :auto);
CSV.write(fname,updf);
