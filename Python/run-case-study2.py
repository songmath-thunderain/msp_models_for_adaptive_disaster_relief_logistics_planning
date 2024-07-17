import os
import time

command = []

instance_option = [" 2 "];

cost_structure = [" 1 ", " 2 "]

nu_levels = [" 1 ", " 2 ", " 5 ", " 10 "];

fc_levels = [" 0.25 ", " 0.5 ", " 2 ", " 4 "]

st_levels = [" 1 ", " 3 ", " 4 "]

for k1 in instance_option:
    for k2 in cost_structure:
        if k2 == " 1 ":
            for k3 in nu_levels:
                command = command + ["python main.py -p solveParams.yaml -d 0 -a 0 " + "-k 1 -o 1000 -ni 3 -nj 10 -t " + k3 + "-s 4 -i " + k1 + "-w 1 -c 1 " + "-st 2 " + "-trans 10 "];
        if k2 == " 2 ":
            command = command + ["python main.py -p solveParams.yaml -d 0 -a 0 " + "-k 1 -o 1000 -ni 3 -nj 10 -t 1000000 " + "-s 4 -i " + k1 + "-w 1 -c 2 " + "-st 2 " + "-trans 10 "];

    for k3 in fc_levels:
        command = command + ["python main.py -p solveParams.yaml -d 0 -a 0 " + "-k 1 -o 1000 -ni 3 -nj 10 -t 1000000 " + "-s 4 -i " + k1 + "-w 1 -c 2 " + "-st 2 " + "-fc " + k3 + "-trans 10 "];
    
    for k3 in st_levels:
        command = command + ["python main.py -p solveParams.yaml -d 0 -a 0 " + "-k 1 -o 1000 -ni 3 -nj 10 -t 1000000 " + "-s 4 -i " + k1 + "-w 1 -c 2 " + "-st " + k3 + "-trans 10 "];

for i in range(len(command)):
    print(command[i]);
    os.system(command[i]);
	
