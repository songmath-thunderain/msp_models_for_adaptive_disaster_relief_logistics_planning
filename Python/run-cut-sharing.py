import os
import time

command = []
ni_options = [" 3 ", " 6 ", " 9 "]
nj_options = [" 10 ", " 20 ", " 30 "]
instance_option = [" -1 ", " 0 "];
tau = [" 0 ", " 0.5 ", " 5 "];
solution_option = [" 0 ", " 1 "];

for k1 in instance_option:
    for ni in ni_options:
        for nj in nj_options:
            for t in tau:
                for k2 in solution_option:
                    command = command + ["python main-cutshare-temp.py -p solveParams.yaml -d 0 -a 0 " + "-k 57 -o 1000 -ni " + ni + "-nj " + nj + " -t " + t + " -i " + k1 + "-w 1 -c 0 -s " + k2];	

for i in range(len(command)):
    print(command[i]);
    os.system(command[i]);
	
