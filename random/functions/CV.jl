#clairvoyance model functions 

#Define the model
function deterministic_model()
 
    #######################
    #Define the model.
    m = Model(optimizer_with_attributes(() -> Gurobi.Optimizer(GRB_ENV), "OutputFlag" => 0, "Presolve" => 0));

    #######################
    #Define the variables.
    @variables(m,
            begin
               0 <= x[i=1:Ni,t=1:T] <= x_cap[i]
               0 <= f[i=1:N0,ii=1:Ni,t=1:T]
               0 <= y[i=1:Ni,j=1:Nj,t=1:T]
               0 <= z[j=1:Nj,t=1:T]
               0 <= v[i=1:Ni,t=1:T]
            end
          );

    #######################
    #Define the objective.
    @objective(m,
               Min,
               sum(sum(sum(cb[i,ii,t]*f[i,ii,t] for ii=1:Ni) for i=1:N0) for t=1:T)
              +sum(sum(ch[i,t]*x[i,t] for i=1:Ni) for t=1:T)
              +sum(sum(f[N0,i,t] for i=1:Ni)*h[t] for t=1:T)
              +sum(sum(sum(ca[i,j,t]*y[i,j,t] for j=1:Nj) for i=1:Ni) for t=1:T)
              +sum(sum(z[j,t] for j=1:Nj)*p for t=1:T)
              +sum(sum(v[i,t] for i=1:Ni)*q for t=1:T)       
               );

    #######################
    #Define the constraints.
    dCons = Dict(); #a dictonary to store all the demand constraint
    for t=1:T
        for i=1:Ni
            if t == 1
                @constraint(m, 
                            x[i,t]
                           +sum(f[i,j,t] for j=1:Ni if j != i)
                           -sum(f[j,i,t] for j=1:N0 if j != i)
                           +sum(y[i,j,t] for j=1:Nj)
                           +v[i,t]
                           == x_0[i]
                            );
               @constraint(m, sum(f[i,j,t] for j=1:Ni if j != i) <= x_0[i]);
            else
                @constraint(m, 
                            x[i,t]
                           +sum(f[i,j,t] for j=1:Ni if j != i)
                           -sum(f[j,i,t] for j=1:N0 if j != i)
                           +sum(y[i,j,t] for j=1:Nj)   
                           +v[i,t]
                           == x[i,t-1]
                            );
               @constraint(m, sum(f[i,j,t] for j=1:Ni if j != i) <= x[i,t-1]);                
            end
        end
        for j=1:Nj
           dCons[t,j] = @constraint(m, z[j,t]+sum(y[i,j,t] for i=1:Ni) >= 0);
        end
        
    end

    return m, x, f, y, z, v, dCons
end

