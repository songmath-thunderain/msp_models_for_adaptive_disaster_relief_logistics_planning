import yaml
import argparse
import sys
from dataClass import *
from CV import *
from TwoStageSP import *
from FA import *

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--solveparam", help="solve parameter yaml file")
    parser.add_argument("-d", "--dissipate_option", type = int, choices = [0,1], help = "dissipate_option")
    parser.add_argument("-a", "--absorbing_option", type = int, choices = [0,1], help = "absorbing_option")
    parser.add_argument("-k", "--k_init", type = int, help = "k_init parameter")
    parser.add_argument("-o", "--oos", type = int, help = "number of out-of-sample scenarios")
    parser.add_argument("-ni", "--Ni", type = int, choices = [3,6,9], help = "number of SPs")
    parser.add_argument("-nj", "--Nj", type = int, choices = [10,20,30], help = "number of DPs")
    parser.add_argument("-t", "--tau", type = float, help = "cost-scaling factor")
    parser.add_argument("-s", "--solution_option", type = int, choices = [0,1,2,3,4], help = "solution options: 0. CV, 1. FA, 2. static2SSP, 3. RH2SSP, 4. WS")
    parser.add_argument("-i", "--instance_option", type = int, choices = [0,1], help = "instance option: 0 -- Synthetic, 1 -- Case Study")
    args = parser.parse_args()
    solveparam_file = args.solveparam
    dissipate_option = args.dissipate_option
    # dissipate_option = 1: has dissipation state, when hurricane intensity reaches state 0, it is absorbed
    # dissipate_option = 0: has no dissipation state, won't be absorbed because of intensity
    absorbing_option = args.absorbing_option
    # absorbing_option = 0: do not allow shipping at the period that hurricane reaches an absorbing state
    # absorbing_option = 1: still allows shipping at the period that hurricane reaches an absorbing state 
    instance_option = args.instance_option
    k_init = args.k_init
    if instance_option == 1:
        # case-study instance, k_init is always 1 regardless
        k_init = 1;
        if dissipate_option == 1:
            print("ERROR! dissipatte_option has to be 0 if running the case study instance!")
            exit(0);
    oos = args.oos
    Ni = args.Ni
    Nj = args.Nj
    tau = args.tau

    with open(solveparam_file, "r") as f:
        params = yaml.safe_load(f)
        max_iter = params['MAX_ITER']
        stall = params['STALL']
        cutviol_maxiter = params['CUTVIOL_MAXITER']
        time_limit = params['TIME_LIMIT']
        cutviol = params['CUTVIOL']

    inputParams = inputParams(dissipate_option, absorbing_option, k_init, oos);
    solveParams = solveParams(max_iter, stall, cutviol_maxiter, time_limit, cutviol);

    hurricaneInstance = hurricaneData();
    networkInstance = networkData(Ni,Nj);

    if instance_option == 0:
        # synthetic instance family
        intensityFile = 'data/synthetic/intensity.csv';
        locationFile = 'data/synthetic/location.csv';
        landfallFile = 'data/synthetic/landfall.csv';
        
        hurricaneInstance.input_from_Syn(intensityFile, locationFile, landfallFile, inputParams)

        netNodesFile = 'data/synthetic/nodes.csv';
        netParamsFile = 'data/synthetic/netParams.csv';
        networkInstance.input_from_Syn(tau,netNodesFile,netParamsFile,hurricaneInstance)

        osfname = "./data/synthetic/OOS" + str(inputParams.k_init) + ".csv"
    elif instance_option == 1:
        # case-study instance
        intensityFile = 'data/case-study/mc_int_transition_prob.csv';
        trackProbFile = 'data/case-study/mc_track_transition_prob_at_t';
        trackErrorFile = 'data/case-study/mc_track_mean_error_at_t';
        landfallFile = 'data/case-study/landfall_deterministic.csv';
        #landfallFile = 'data/case-study/landfall.csv';
        hurricaneInstance.input_from_Case(intensityFile, trackProbFile, trackErrorFile, landfallFile);

        netFolderPath = 'data/case-study';
        netParamsFile = 'data/case-study/netParams.csv';
        networkInstance.input_from_Case(tau, netFolderPath, netParamsFile, hurricaneInstance);

        osfname = "./data/case-study/OOS1_deterministic.csv"
        #osfname = "./data/case-study/OOS1.csv"
    else:
        print("ERROR: instance_option has to be 0 or 1!")
        exit(0);

    option = args.solution_option
    if option == 0:
        CV = CV(inputParams,hurricaneInstance,networkInstance)
        CV.clairvoyant_eval(osfname) 
    elif option == 1:
        FA = FA(inputParams,solveParams,hurricaneInstance,networkInstance)
        FA.FOSDDP_eval(osfname)
    elif option == 2:
        TwoStageSP = TwoStageSP(inputParams,solveParams,hurricaneInstance,networkInstance)
        TwoStageSP.static_2SSP_eval(osfname)
    elif option == 3:
        TwoStageSP = TwoStageSP(inputParams,solveParams,hurricaneInstance,networkInstance)
        TwoStageSP.RH_2SSP_eval(osfname)
    elif option == 4:
        TwoStageSP = TwoStageSP(inputParams,solveParams,hurricaneInstance,networkInstance)
        TwoStageSP.WS_eval(osfname)
    else:
        print("This option is not available!")
        sys.exit(0);