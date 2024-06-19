import os
import time

command = []

instance_option = [" 1 ", " 2 "];
cost_structure = [" 0 ", " 1 ", " 2 "];
tau = [" 2 ", " 5 ", " 10 "]

for k1 in instance_option:
    for k2 in cost_structure:
        if k2 == " 0 ":
            command = command + ["python main.py -p solveParams.yaml -d 0 -a 0 " + "-k 1 -o 1000 -ni 3 -nj 10 -t 0 " + "-s 5 -i " + k1 + "-w 1 -c " + k2];
        if k2 == " 1 ":
            for k3 in tau:
                command = command + ["python main.py -p solveParams.yaml -d 0 -a 0 " + "-k 1 -o 1000 -ni 3 -nj 10 -t " + k3 + "-s 5 -i " + k1 + "-w 1 -c " + k2 + "-st 2 "];
        if k2 == " 2 ":
            command = command + ["python main.py -p solveParams.yaml -d 0 -a 0 " + "-k 1 -o 1000 -ni 3 -nj 10 -t 1000000 " + "-s 5 -i " + k1 + "-w 1 -c " + k2 + "-st 2 "];

for i in range(len(command)):
    print(command[i]);
    os.system(command[i]);
	
