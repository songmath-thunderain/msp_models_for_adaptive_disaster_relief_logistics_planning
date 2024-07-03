import os
import time

command = []
ni_options = [" 3 ", " 6 " , " 9 "]
nj_options = [" 10 ", " 20 ", " 30 "]
instance_option = [" -1 ", " 0 "];
tau2 = [" 2 ", " 5 ", " 10 "];
tau = [" 0 ", " 0.5 ", " 5 "]
for k1 in instance_option:
    for ni in ni_options:
        for nj in nj_options:
            for t in tau:
                command = command + ["python main.py -p solveParams.yaml -d 0 -a 0 " + "-k 57 -o 1000 -ni " + ni + "-nj " + nj + " -t " + t + "-s 1 -i " + k1 + "-w 1 -c 0 "];	
    for ni in ni_options:
        for nj in nj_options:
            for t in tau2:
                command = command + ["python main.py -p solveParams.yaml -d 0 -a 0 " + "-k 57 -o 1000 -ni " + ni + "-nj " + nj + " -t " + t + "-s 1 -i " + k1 + "-w 1 -c 1 " + "-st 2 "];	

for i in range(len(command)):
    print(command[i]);
    os.system(command[i]);
	
