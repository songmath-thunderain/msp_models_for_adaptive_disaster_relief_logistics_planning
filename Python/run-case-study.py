import os
import time

command = []

absorbing_option = [" 0 ", " 1 "];
instance_option = [" 1 ", " 2 "];
cost_structure = [" 0 ", " 1 "];
tau = [" 0 ", " 0.5 ", " 5 "];
tau2 = [" 2 ", " 5 ", " 10 "];

for k1 in instance_option:
    for k2 in cost_structure:
        for k3 in absorbing_option:
            if k2 == " 0 ":
                for t in tau:
                    command = command + ["python main.py -p solveParams.yaml -d 0 -a " + k3 + "-k 1 -o 1000 -ni 3 -nj 10 -t " + t + "-s 5 -i " + k1 + "-w 1 -c " + k2];	
            if k2 == " 1 ":
                for t in tau2:
                    command = command + ["python main.py -p solveParams.yaml -d 0 -a " + k3 + "-k 1 -o 1000 -ni 3 -nj 10 -t " + t + "-s 5 -i " + k1 + "-w 1 -c " + k2 + "-st 2"];	

for i in range(len(command)):
    print(command[i]);
    os.system(command[i]);
	
