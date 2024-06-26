
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
objs_fa, fa_bar, fa_low, fa_high, elapsed_fa = FOSDDP_eval_offline();

fname = "./output/benchmark/FAresults.csv"
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
