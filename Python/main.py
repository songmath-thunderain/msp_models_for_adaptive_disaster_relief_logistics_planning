import yaml
import argparse
import sys
import csv
import time
import pickle
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
    parser.add_argument("-c", "--cost_structure", type = int, choices = [-1,0,1], help = "cost structure option: -1. safe time but time dependent only, 0. time increasing, 1. safe time but state dependent")
    parser.add_argument("-t", "--tau", type = float, help = "cost-scaling factor")
    parser.add_argument("-st", "--safe_time", required = False, type = int, help = "safe time parameter to determine cost surge")
    parser.add_argument("-s", "--solution_option", type = int, choices = [0,1,2,3,4,5], help = "solution options: 0. CV, 1. FA, 2. static2SSP, 3. RH2SSP, 4. WS, 5. All")
    parser.add_argument("-i", "--instance_option", type = int, choices = [-1,0,1], help = "instance option: -1. Synthetic D-landfall, 0. Synthetic R-landfall, 1. Case Study D-landfall")
    parser.add_argument("-w", "--write_option", type = int, choices = [0,1], help = "0. do not write to CSV, 1. write results to CSV")
    args = parser.parse_args()
    solveparam_file = args.solveparam
    dissipate_option = args.dissipate_option
    # dissipate_option = 1: has dissipation state, when hurricane intensity reaches state 0, it is absorbed
    # dissipate_option = 0: has no dissipation state, won't be absorbed because of intensity
    absorbing_option = args.absorbing_option
    # absorbing_option = 0: do not allow shipping at the period that hurricane reaches an absorbing state
    # absorbing_option = 1: still allows shipping at the period that hurricane reaches an absorbing state 
    instance_option = args.instance_option
    write_option = args.write_option
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
    cost_structure = args.cost_structure
    # cost_structure = -1: modified cost structure but in the deterministic landfall case
    # cost_structure = 0: original cost structure: cost increases by tau in every period
    # cost_structure = 1: modified cost structure: cost stays the same until the expected landfall time is within the safe time threshold
    # safe_time is optional
    safe_time = None
    if args.safe_time is not None:
        safe_time = args.safe_time
    else:
        if cost_structure == 1 or cost_structure == -1:
            print("Error! safe_time parameter is not defined for cost_structure == -1 or 1!")
            exit(0);

    with open(solveparam_file, "r") as f:
        params = yaml.safe_load(f)
        max_iter = params['MAX_ITER']
        stall = params['STALL']
        cutviol_maxiter = params['CUTVIOL_MAXITER']
        time_limit = params['TIME_LIMIT']
        cutviol = params['CUTVIOL']

    inputParams = inputParams(dissipate_option, absorbing_option, cost_structure, safe_time, tau, k_init, oos);
    solveParams = solveParams(max_iter, stall, cutviol_maxiter, time_limit, cutviol);

    hurricaneInstance = hurricaneData();
    networkInstance = networkData(Ni,Nj);

    print("Reading data...");
    start_time = time.time()

    ISpaths = None
    if instance_option == -1:
        # synthetic instance family: deterministic landfall time
        intensityFile = 'data/synthetic/intensity.csv';
        locationFile = 'data/synthetic/location.csv';
        landfallFile = 'data/synthetic/landfall-D.csv';
        
        hurricaneInstance.input_from_Syn(intensityFile, locationFile, landfallFile, inputParams)

        netNodesFile = 'data/synthetic/nodes.csv';
        netParamsFile = 'data/synthetic/netParams.csv';
        networkInstance.input_from_Syn(cost_structure,safe_time,tau,netNodesFile,netParamsFile,hurricaneInstance)

        osfname = "./data/synthetic/OOS" + str(inputParams.k_init) + "-D.csv"

    elif instance_option == 0:
        # synthetic instance family: random landfall time
        intensityFile = 'data/synthetic/intensity.csv';
        locationFile = 'data/synthetic/location.csv';
        landfallFile = 'data/synthetic/landfall.csv';
        
        hurricaneInstance.input_from_Syn(intensityFile, locationFile, landfallFile, inputParams)

        netNodesFile = 'data/synthetic/nodes.csv';
        netParamsFile = 'data/synthetic/netParams.csv';
        networkInstance.input_from_Syn(cost_structure,safe_time,tau,netNodesFile,netParamsFile,hurricaneInstance)

        osfname = "./data/synthetic/OOS" + str(inputParams.k_init) + ".csv"

        if cost_structure == 1:
            ISfile = open('data/synthetic/in_sample_100.dat', 'rb')
            ISpaths = pickle.load(ISfile)
            ISfile.close()

    elif instance_option == 1:
        # case-study instance: deterministic landfall time
        '''
        # old case study input interface
        intensityFile = 'data/case-study/mc_int_transition_prob.csv';
        trackProbFile = 'data/case-study/mc_track_transition_prob_at_t';
        trackErrorFile = 'data/case-study/mc_track_mean_error_at_t';
        landfallFile = 'data/case-study/landfall_deterministic.csv';
        #landfallFile = 'data/case-study/landfall.csv';
        hurricaneInstance.input_from_Case(intensityFile, trackProbFile, trackErrorFile, landfallFile);

        netFolderPath = 'data/case-study';
        netParamsFile = 'data/case-study/netParams.csv';
        networkInstance.input_from_Case(cost_structure,safe_time,tau,netFolderPath,netParamsFile,hurricaneInstance);

        osfname = "./data/case-study/OOS1_deterministic.csv"
        #osfname = "./data/case-study/OOS1.csv"

        if cost_structure == 1:
            ISfile = open('data/case-study/in_sample_100.dat', 'rb')
            ISpaths = pickle.load(ISfile)
            ISfile.close()
        '''
        # new case study input interface

        forecastFile = 'data/case-study/SC-network/Florence_forecast.csv';
        MCFile = 'data/case-study/SC-network/deterministic/pi_mssp_d.json';
        hurricaneInstance.input_from_Case_new(forecastFile, MCFile);

        netFolderPath = 'data/case-study/SC-network';
        netParamsFile = 'data/case-study/SC-network/netParams.csv';
        networkInstance.input_from_Case_new(cost_structure,safe_time,tau,netFolderPath,netParamsFile,hurricaneInstance);

    else:
        print("ERROR: instance_option has to be -1, 0 or 1!")
        exit(0);

    elapse_time = time.time() - start_time
    print("Elapse time = ", elapse_time);

    option = args.solution_option
    if option == 0 or option == 5:
        CV = CV(inputParams,hurricaneInstance,networkInstance)
        obj, CI, elapsed_time = CV.clairvoyant_eval(osfname)
        if write_option == 1:
            with open('output/CVresults.csv', 'a') as myfile:
                writer = csv.writer(myfile, delimiter =',')
                writer.writerow([instance_option,dissipate_option,absorbing_option,k_init,Ni,Nj,tau,obj,CI,elapsed_time])
    if option == 1 or option == 5:
        FA = FA(inputParams,solveParams,hurricaneInstance,networkInstance)
        obj, CI, train_time, test_time = FA.FOSDDP_eval(osfname)
        if write_option == 1:
            with open('output/FAresults.csv', 'a') as myfile:
                writer = csv.writer(myfile, delimiter =',')
                writer.writerow([instance_option,dissipate_option,absorbing_option,k_init,Ni,Nj,tau,obj, CI, train_time, test_time])
    if option == 2 or option == 5:
        TwoStageSP = TwoStageSP(inputParams,solveParams,hurricaneInstance,networkInstance,ISpaths)
        obj, CI, train_time, test_time = TwoStageSP.static_2SSP_eval(osfname)
        if write_option == 1:
            with open('output/static2SSPresults.csv', 'a') as myfile:
                writer = csv.writer(myfile, delimiter =',')
                writer.writerow([instance_option,dissipate_option,absorbing_option,k_init,Ni,Nj,tau,obj, CI, train_time, test_time])
    if option == 3 or option == 5:
        if option == 3:
            TwoStageSP = TwoStageSP(inputParams,solveParams,hurricaneInstance,networkInstance,ISpaths)
        obj, CI, elapsed_time = TwoStageSP.RH_2SSP_eval(osfname)
        if write_option == 1:
            with open('output/rolling2SSPresults.csv', 'a') as myfile:
                writer = csv.writer(myfile, delimiter =',')
                writer.writerow([instance_option,dissipate_option,absorbing_option,k_init,Ni,Nj,tau,obj,CI,elapsed_time])
    if option == 4 or option == 5:
        if option == 4:
            TwoStageSP = TwoStageSP(inputParams,solveParams,hurricaneInstance,networkInstance,ISpaths)
        obj, CI, train_time, test_time = TwoStageSP.WS_eval(osfname)
        if write_option == 1:
            with open('output/WSresults.csv', 'a') as myfile:
                writer = csv.writer(myfile, delimiter =',')
                writer.writerow([instance_option,dissipate_option,absorbing_option,k_init,Ni,Nj,tau,obj, CI, train_time, test_time])