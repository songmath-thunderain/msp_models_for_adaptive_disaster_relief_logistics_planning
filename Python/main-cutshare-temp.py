import yaml
import argparse
import sys
import csv
import time
import pickle
from dataClass import *
from CV import *
from TwoStageSP import *
from FACutshareTemp import *

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--solveparam", help="solve parameter yaml file")
    parser.add_argument("-d", "--dissipate_option", type = int, choices = [0,1], help = "dissipate_option")
    parser.add_argument("-a", "--absorbing_option", type = int, choices = [0,1], help = "absorbing_option")
    parser.add_argument("-k", "--k_init", type = int, help = "k_init parameter")
    parser.add_argument("-o", "--oos", type = int, help = "number of out-of-sample scenarios")
    parser.add_argument("-ni", "--Ni", type = int, choices = [3,6,9], help = "number of SPs")
    parser.add_argument("-nj", "--Nj", type = int, choices = [10,20,30], help = "number of DPs")
    parser.add_argument("-c", "--cost_structure", type = int, choices = [0,1,2], help = "cost structure option: 0. time increasing, 1. safe time with logistics cost surge, 2. only price surge by tau in logistics cost")
    parser.add_argument("-t", "--tau", type = float, help = "cost-scaling factor")
    parser.add_argument("-st", "--safe_time", required = False, type = int, help = "safe time parameter to determine cost surge")
    parser.add_argument("-s", "--solution_option", type = int, choices = [0,1], help = "solution options: 0. no cut sharing, 1. has cut sharing")
    parser.add_argument("-i", "--instance_option", type = int, choices = [-1,0,1,2], help = "instance option: -1. Synthetic D-landfall, 0. Synthetic R-landfall, 1. Case Study D-landfall, 2. Case Study R-landfall")
    parser.add_argument("-w", "--write_option", type = int, choices = [0,1], help = "0. do not write to CSV, 1. write results to CSV")
    parser.add_argument("-fc", "--flow_capacity", type = float, choices = [0.125,0.25,0.5,1,2,4], required = False, help = "flow capacity level: needs to be 0.125, 0.25, 0.5, 1, 2, or 4")
    parser.add_argument("-pen", "--penalty", type = float, choices = [0.25,0.5,1,2,4], required = False, help = "penalty level: needs to be 0.25, 0.5, 1, 2, or 4")
    parser.add_argument("-inv", "--inventory", type = float, choices = [0,0.5,1,2], required = False, help = "inventory level: needs to be 0, 0.5, 1, or 2")
    parser.add_argument("-trans", "--transportation", type = float, choices = [1,5,10,20], required = False, help = "transportation level: needs to be 1,5,10,or 20")

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
    if instance_option == 1 or instance_option == 2:
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
    # cost_structure = 0: original cost structure: cost increases by tau in every period
    # cost_structure = 1: modified cost structure: cost stays the same until the expected landfall time is within the safe time threshold, then the price surge by tau
    # cost_structure = 2: only price surge by tau in logistics cost 
    # safe_time is optional
    safe_time = None
    if args.safe_time is not None:
        safe_time = args.safe_time
    else:
        if cost_structure == 1 or cost_structure == 2:
            print("Error! safe_time parameter is not defined for cost_structure == 1 or 2!")
            exit(0);
        if args.solution_option == -1:
            print("Error! safe_time parameter is not defined for solution_option == -1, i.e., naiveWS!")
            exit(0);
    arc_option = False
    if instance_option == 1 or instance_option == 2:
        # case study, always assume the existence of arc.xlsx file
        arc_option = True
    # arc_option = 0: no arc restriction is imposed
    # arc_option = 1: arc restriction is imposed according to a file that specifies the set of allowed arcs (currently only between SPs and DPs)
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

    fc_level = 1;
    if args.flow_capacity != None:
        fc_level = args.flow_capacity;
    
    penalty_level = 1;
    if args.penalty != None:
        penalty_level = args.penalty;
    
    inventory_level = 1;
    if args.inventory != None:
        inventory_level = args.inventory;
    
    transportation_level = 1;
    if args.transportation != None:
        transportation_level = args.transportation;

    ISpaths = None
    if instance_option == -1:
        # synthetic instance family: deterministic landfall time
        intensityFile = 'data/synthetic/intensity.csv';
        locationFile = 'data/synthetic/location.csv';
        landfallFile = 'data/synthetic/landfall-D.csv';
        
        hurricaneInstance.input_from_Syn(intensityFile, locationFile, landfallFile, inputParams)

        netNodesFile = 'data/synthetic/nodes.csv';
        netParamsFile = 'data/synthetic/netParams.csv';
        networkInstance.input_from_Syn(cost_structure,safe_time,tau,netNodesFile,netParamsFile,hurricaneInstance,arc_option,0)

        osfname = "./data/synthetic/OOS" + str(inputParams.k_init) + "-D.csv"

    elif instance_option == 0:
        # synthetic instance family: random landfall time
        intensityFile = 'data/synthetic/intensity.csv';
        locationFile = 'data/synthetic/location.csv';
        landfallFile = 'data/synthetic/landfall.csv';
        
        hurricaneInstance.input_from_Syn(intensityFile, locationFile, landfallFile, inputParams)

        netNodesFile = 'data/synthetic/nodes.csv';
        netParamsFile = 'data/synthetic/netParams.csv';
        networkInstance.input_from_Syn(cost_structure,safe_time,tau,netNodesFile,netParamsFile,hurricaneInstance,arc_option,1)

        osfname = "./data/synthetic/OOS" + str(inputParams.k_init) + ".csv"

        if cost_structure == 1 or cost_structure == 2:
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
        absorbingFile = None; # deterministic landfall case: no need to supply the absorbingFile
        MCFile = 'data/case-study/SC-network/deterministic/pi_mssp_d.json';
        hurricaneInstance.input_from_Case_new(absorbingFile, MCFile);

        netFolderPath = 'data/case-study/SC-network/';
        networkInstance.input_from_Case_new(cost_structure,safe_time,tau,netFolderPath,hurricaneInstance,arc_option,0,fc_level,penalty_level,inventory_level,transportation_level);
 
        osfname = "./data/case-study/SC-network/deterministic/OOS" + str(inputParams.k_init) + "-D.csv"

    elif instance_option == 2:
        # new case study input interface
        absorbingFile = 'data/case-study/SC-network/random/absorb_mssp_r_fix_along_True.json'; 
        MCFile = 'data/case-study/SC-network/random/pi_mssp_r_fix_along_True.json';
        hurricaneInstance.input_from_Case_new(absorbingFile, MCFile);

        print("states = ", hurricaneInstance.states);
        print("absorbing states = ", hurricaneInstance.absorbing_states);

        netFolderPath = 'data/case-study/SC-network/';
        networkInstance.input_from_Case_new(cost_structure,safe_time,tau,netFolderPath,hurricaneInstance,arc_option,1,fc_level,penalty_level,inventory_level,transportation_level);

        osfname = "./data/case-study/SC-network/random/OOS" + str(inputParams.k_init) + ".csv"

        if cost_structure == 1 or cost_structure == 2:
            ISfile = open('data/case-study/SC-network/random/in_sample_100.dat', 'rb')
            ISpaths = pickle.load(ISfile)
            ISfile.close()

    else:
        print("ERROR: instance_option has to be -1, 0, 1, or 2!")
        exit(0);

    elapse_time = time.time() - start_time
    print("Elapse time = ", elapse_time);

    outputpath = None;
    if instance_option == -1 or instance_option == 0:
        # synthetic instance, output path is set to synthetic
        outputpath = 'output/synthetic/';
    elif instance_option == 1 or instance_option == 2:
        # case study instance, output path is set to case study
        outputpath = 'output/case-study/';
    else:
        print("ERROR: instance_option has to be -1, 0, 1, or 2!") 
        exit(0);

    option = args.solution_option
    if safe_time is None:
        safe_time = 0; # just print out something trivial
    
    FA = FA(inputParams,solveParams,hurricaneInstance,networkInstance,option)
    [LB, iter, obj, CI, train_time, test_time] = FA.FOSDDP_eval(osfname)
    if write_option == 1:
        with open(outputpath+'FAresults-cutsharing.csv', 'a') as myfile:
            writer = csv.writer(myfile, delimiter =',')
            writer.writerow([instance_option,cost_structure,dissipate_option,absorbing_option,k_init,Ni,Nj,tau,safe_time,option,LB,iter,obj,CI,train_time,test_time])
