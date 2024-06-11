import os
import time

command = []
ni_options = [" 3 ", " 6 " , " 9 "]
nj_options = [" 10 ", " 20 ", " 30 "]
instance_option = [" -1 ", " 0 "];
cost_structure = [" 0 ", " 1 "];
tau = [" 0 ", " 0.5 ", " 5 "];
tau2 = [" 2 ", " 5 ", " 10 "];

for k1 in instance_option:
    for ni in ni_options:
        for nj in nj_options:
            for k2 in cost_structure:
                    if k2 == " 0 ":
                        for t in tau:
                            command = command + ["python main.py -p solveParams.yaml -d 0 -a 0 " + "-k 57 -o 1000 -ni " + ni + "-nj " + nj + " -t " + t + "-s 5 -i " + k1 + "-w 1 -c " + k2];	
                    if k2 == " 1 ":
                        for t in tau2:
                            command = command + ["python main.py -p solveParams.yaml -d 0 -a 0 " + "-k 57 -o 1000 -ni " + ni + "-nj " + nj + " -t " + t + "-s 5 -i " + k1 + "-w 1 -c " + k2 + "-st 2 "];	

for i in range(len(command)):
    print(command[i]);
    os.system(command[i]);
	
