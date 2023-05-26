#clairvoyant solution 

start=time();
osfname = "./data/OOS"*string(k_init)*".csv";
OS_paths = Matrix(CSV.read(osfname,DataFrame)); #read the out-of-sample file
# we do not have this second layer now [REVISION]
# OS_M = Matrix(CSV.read("./data/inOOS.csv",DataFrame))[:,1] #read the second layer OOS
objs_cv = zeros(nbOS);

xval_cv = Array{Any,1}(undef,nbOS); 
fval_cv = Array{Any,1}(undef,nbOS);
yval_cv = Array{Any,1}(undef,nbOS); 
zval_cv = Array{Any,1}(undef,nbOS); 
vval_cv = Array{Any,1}(undef,nbOS);

for s=1:nbOS
	#find the period when the hurricane made landfall && intensity != 1
	#println("OS_paths[", s, "] = ", OS_paths[s,1:T]);
	absorbingT = -1;
	if dissipate_option == 1
		absorbingT = findfirst(x -> (S[x][3] == Nc || S[x][1] == 1), OS_paths[s,1:T]);
	else
		absorbingT = findfirst(x -> (S[x][3] == Nc), OS_paths[s,1:T]);
	end
	if S[OS_paths[s,absorbingT]][1] == 1
		# we knew nothing is going to happen, so cost would be 0
		continue;
	else
		#define the clairvoyant model 
		m_cv, x_cv, f_cv, y_cv, z_cv, v_cv, dCons_cv = deterministic_model();
		k_t = OS_paths[s,absorbingT] # state corresponding to the landfall
		# m = OS_M[s]; [REVISION]
		for j=1:Nj, t=1:T
			if t == absorbingT
				set_normalized_rhs(dCons_cv[absorbingT,j], SCEN[k_t][j]); #set the RHS of demand constraint
				if absorbing_option == 0
					for i=1:N0
						for ii=1:Ni
							set_upper_bound(f_cv[i,ii,t],0);
						end
					end
				end
			else
				set_normalized_rhs(dCons_cv[t,j], 0); #set the RHS of demand constraint
			end
		end
		 
		optimize!(m_cv); #solve the model        
		status = termination_status(m_cv); #check the status 
		if status != MOI.OPTIMAL
			println(" Model in clairvoyant is ", status);
			exit(0);
		else
			xval_cv[s] = value.(x_cv);
			fval_cv[s] = value.(f_cv);                 
			yval_cv[s] = value.(y_cv);
			zval_cv[s] = value.(z_cv);
			vval_cv[s] = value.(v_cv);
			objs_cv[s] = objective_value(m_cv);
		end
	end
end


cv_bar = mean(objs_cv);
cv_std = std(objs_cv);
cv_low = cv_bar-1.96*cv_std/sqrt(nbOS);
cv_high = cv_bar+1.96*cv_std/sqrt(nbOS);
println("Clairvoyant....");
println("μ ± 1.96*σ/√NS = ", cv_bar, "±", [cv_low,cv_high]);
elapsed_cv = time() - start;
vals = [xval_cv, fval_cv, yval_cv, zval_cv, vval_cv];

fname = "./output/benchmark/CVresults.csv"
df = CSV.read(fname,DataFrame);
results_cv = Matrix(df);
results_cv[inst,1] = 0
results_cv[inst,2] = cv_bar
results_cv[inst,3] = cv_bar-cv_low
results_cv[inst,4] = cv_bar+cv_low
results_cv[inst,5] = elapsed_cv
results_cv[inst,6] = 0

updf = DataFrame(results_cv, :auto);
CSV.write(fname,updf)
