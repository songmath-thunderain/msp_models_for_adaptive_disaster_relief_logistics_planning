import yaml
import argparse
import sys
from dataInput import *
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
    parser.add_argument("-s", "--solution_option", type = int, choices = [1,2,3,4], help = "solution options: 1. FA, 2. static2SSP, 3. RH2SSP, 4. WS")
    args = parser.parse_args()
    solveparam_file = args.solveparam
    dissipate_option = args.dissipate_option
    absorbing_option = args.absorbing_option
    k_init = args.k_init
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

    intensityFile = 'data/synthetic/intensity.csv';
    locationFile = 'data/synthetic/location.csv';
    landfallFile = 'data/synthetic/landfall_7.csv';

    hurricaneDataSet = hurricaneInput(intensityFile, locationFile, landfallFile, inputParams)

    netNodesFile = 'data/synthetic/nodes.csv';
    netParamsFile = 'data/synthetic/netParams.csv';
    networkDataSet = networkInput(Ni,Nj,tau,netNodesFile,netParamsFile,hurricaneDataSet)

    osfname = "./data/synthetic/OOS" + str(inputParams.k_init) + ".csv"

    option = args.solution_option
    if option == 1:
        FA = FA(inputParams,solveParams,hurricaneDataSet,networkDataSet)
        FA.FOSDDP_eval(osfname)
    elif option == 2:
        static_2SSP_eval(networkDataSet, hurricaneDataSet, inputParams, solveParams, osfname)
    elif option == 3:
        RH_2SSP_eval(networkDataSet, hurricaneDataSet, inputParams, solveParams, osfname)
    elif option == 4:
        WS_eval(networkDataSet, hurricaneDataSet, inputParams, solveParams, osfname)
    else:
        print("This option is not available!")
        sys.exit(0);