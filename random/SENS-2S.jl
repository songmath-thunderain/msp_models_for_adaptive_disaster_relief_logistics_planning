#max_iter = 100;
s = 1;
t_roll = 1;
x_init = x_0;

osfname = "./data/OOS"*string(k_init)*".csv";
OS_paths = Matrix(CSV.read(osfname,DataFrame)); #read the out-of-sample file
#OS_M = Matrix(CSV.read("./data/inOOS.csv",DataFrame))[:,1] #read the second layer OOS [REVISION: no need any more]

#define the the model.
master, x, f, θ, subproblem, y2, xCons, dCons, rCons = RH_2SSP_define_models(t_roll,OS_paths[s,t_roll],x_init);

#solve the model.
LB_st2SSP, UB_st2SSP, xval_st2SSP, fval_st2SSP, θval_st2SSP = RH_2SSP_solve_roll(OS_paths[s,t_roll],t_roll,master,subproblem,x,f,θ,y2,xCons,dCons,rCons);

procurmnt_amount = zeros(T); 

for t = 1:T
	procurmnt_amount[t] += sum(fval_st2SSP[N0,i,t] for i=1:Ni)
end

fname = "./output/2S-sensitivity.csv"
df = CSV.read(fname,DataFrame);
results_fa = Matrix(df);
results_fa[inst,1:T] = procurmnt_amount
results_fa[inst,end] = inst; 
updf = DataFrame(results_fa, :auto);
CSV.write(fname,updf)
