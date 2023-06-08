#change the currnet working directory to where the file is executed from
cwd = "/"*relpath((@__FILE__)*"/..","/");
cd(cwd) 
#comment this line if don't have bots setup
#include(homedir()*"/misc/tgtext.jl");

#1: instance number 
#2: cost scaling factor
#3: number of scenarios in the out-of-sample 
#4: absorbing option
PARAMS = ARGS;    
instname = PARAMS[1]*"-"*PARAMS[2]*"SPs-"*PARAMS[3]*"DPs-"*PARAMS[4]*"ρ";
#tg_sendtext("Julia: Now starting instance $instname"); #comment this line if don't have bots setup
inst = parse(Int, PARAMS[1]); #instance number 
factor = parse(Float64, PARAMS[2]); #cost scaling factor
nbOS = parse(Int, PARAMS[3]); #number of scenarios in the out-of-sample
absorbing_option = parse(Int, PARAMS[4]); # whether or not we allow MDC/SP operation in the absorbing state
dissipate_option = 0; # We do not treat intensity = 0 as an absorbing state in the case study
FA_option = 1; # We turn on cut sharing

penaltyFactor = parse(Float64, PARAMS[5]); # penalty parameter: 10, 50, 100

k_init = 1;

include("packages.jl");
include("./case-study/case-study.jl");

p = penaltyFactor*base;

include("./functions/functions.jl");


nodeLists = createNodes(k_init); # List of MC states in each stage

nodeScenList, nodeScenWeights = createNodeScens(k_init, nodeLists);

#create_OSpaths(k_init)

#create gurobi environment
const GRB_ENV = Gurobi.Env();

#Clairvoyance solution 
#include("./policy/hurricane_CV.jl");

#Fully adaptive model
#include("./policy/hurricane_FA.jl");
 
#Rolling-horizon 2SSP
#include("./policy/hurricane_rolling_twostage.jl");

#Wait-and-see
#include("./policy/wait_and_see.jl");

#Static 2SSP
#include("./policy/hurricane_static_twostage.jl");

#Deterministic FA 
#include("./policy/hurricane_deterministicFA.jl");

#sensitivity analysis
include("SENS-FA.jl");


#tg_sendtext("Julia: $instname is DONE!"); #comment this line if don't have bots setup
