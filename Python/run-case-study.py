import os
import time

command = []

instance_option = [" 1 ", " 2 "];

cost_factor = [" 1 ", " 2 ", " 5 ", " 10 "];

for k1 in instance_option:
    for k3 in cost_factor:
        command = command + ["python main.py -p solveParams.yaml -d 0 -a 0 " + "-k 1 -o 1000 -ni 3 -nj 10 -t " + k3 + "-s 3 -i " + k1 + "-w 1 -c 1 " + "-st 2 " + "-trans 10 "];
    command = command + ["python main.py -p solveParams.yaml -d 0 -a 0 " + "-k 1 -o 1000 -ni 3 -nj 10 -t 1000000 " + "-s 3 -i " + k1 + "-w 1 -c 2 " + "-st 2 " + "-trans 10 "];

for i in range(len(command)):
    print(command[i]);
    os.system(command[i]);
	
