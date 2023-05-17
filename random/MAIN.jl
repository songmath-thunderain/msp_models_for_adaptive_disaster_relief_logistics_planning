#change the currnet working directory to where the file is executed from
cwd = "/"*relpath((@__FILE__)*"/..","/");
cd(cwd) 
#comment this line if don't have bots setup
#include(homedir()*"/misc/tgtext.jl");

#1: instance number 
#2: number of supply points
#3: number of demand points
#4: cost scaling factor
#5: number of scenarios in the out-of-sample 
#6: initial state
#7: absorbing_option
#8: dissipate_option
PARAMS = ARGS;    
instname = PARAMS[1]*"-"*PARAMS[2]*"SPs-"*PARAMS[3]*"DPs-"*PARAMS[4]*"ρ";
#tg_sendtext("Julia: Now starting instance $instname"); #comment this line if don't have bots setup
inst = parse(Int, PARAMS[1]); #instance number 
Ni = parse(Int, PARAMS[2]); #number of supply points 
Nj = parse(Int, PARAMS[3]); #number of demand points
factor = parse(Float64, PARAMS[4]); #cost scaling factor
nbOS = parse(Int, PARAMS[5]); #number of scenarios in the out-of-sample
k_init = parse(Int, PARAMS[6]); #initial state
absorbing_option = parse(Int, PARAMS[7]); # whether or not we allow MDC/SP operation in the absorbing state
dissipate_option = parse(Int, PARAMS[8]); # whether or not we treat intensity = 0 as an absorbing state

include("packages.jl");
include("./data/data.jl");
include("./functions/functions.jl");

println("Na = ", Na);
println("Nb = ", Nb);
println("Nc = ", Nc);

nodeLists = createNodes(k_init); # List of MC states in each stage
nodeScenList, nodeScenWeights = createNodeScens(k_init, nodeLists);

#create_OSpaths(k_init)

#create gurobi environment
const GRB_ENV = Gurobi.Env();

#Clairvoyance solution 
include("./policy/hurricane_CV.jl");

#Fully adaptive model
include("./policy/hurricane_FA.jl");
 
#Rolling-horizon 2SSP
include("./policy/hurricane_rolling_twostage.jl");

#Wait-and-see
include("./policy/wait_and_see.jl");

#Static 2SSP
include("./policy/hurricane_static_twostage.jl");

#Deterministic FA 
#include("./policy/hurricane_deterministicFA.jl");

#sensitivity analysis
#include("SENS.jl");


#tg_sendtext("Julia: $instname is DONE!"); #comment this line if don't have bots setup
