import os
import time

command = []

instance_option = [" 1 ", " 2 "];
penalty_level = [" 0.125 "]

for k1 in instance_option:
    for k2 in penalty_level:
        command = command + ["python main.py -p solveParams.yaml -d 0 -a 0 " + "-k 1 -o 1000 -ni 3 -nj 10 -t 1000000 " + "-s 5 -i " + k1 + "-w 1 -c 2 " + "-st 2 " + "-fc " + k2];

for i in range(len(command)):
    print(command[i]);
    os.system(command[i]);
	
