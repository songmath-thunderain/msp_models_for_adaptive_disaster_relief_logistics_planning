
#Fully adaptive with deterministic landfall 
#Defne the model 
m_fa, x_fa, f_fa, ϴ_fa, FB1Cons_fa, FB2Cons_fa, model_final, FB_final, dCons_final, y_final = define_models_FAD();


#train the model
LB, train_time, iter = train_models_offline_FAD();

println("***********************************************")
println("***********************************************")
println("finished training!")
println("LB = ", LB[end])
println("training time = ", train_time)
println("number of iterations = ", iter)


#evaluate the model 
objs_fa, fa_bar, fa_low, fa_high, elapsed_fa = FOSDDP_eval_offline_FAD();

fname = "./output/benchmark/deterministicFAresults.csv"
df = CSV.read(fname,DataFrame);
results_fa = Matrix(df);

results_fa[inst,1] = LB[end]

results_fa[inst,2] = fa_bar

results_fa[inst,3] = fa_bar-fa_low

results_fa[inst,4] = train_time

results_fa[inst,5] = elapsed_fa

results_fa[inst,6] = iter

updf = DataFrame(results_fa, :auto);
CSV.write(fname,updf)
