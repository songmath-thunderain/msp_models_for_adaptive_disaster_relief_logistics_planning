import os
import time

command = []

instance_option = [" 1 ", " 2 "];

fc_levels = [" 0.25 ", " 1 ", " 4 "]

st_levels = [" 1 ", " 3 ", " 4 "]

for k1 in instance_option:

    for k3 in fc_levels:
        command = command + ["python main.py -p solveParams.yaml -d 0 -a 0 " + "-k 1 -o 5000 -ni 3 -nj 10 -t 1000000 " + "-s 4 -i " + k1 + "-w 2 -c 2 " + "-st 2 " + "-fc " + k3 + "-trans 10 "];
    
    for k3 in st_levels:
        command = command + ["python main.py -p solveParams.yaml -d 0 -a 0 " + "-k 1 -o 5000 -ni 3 -nj 10 -t 1000000 " + "-s 4 -i " + k1 + "-w 2 -c 2 " + "-st " + k3 + "-trans 10 "];

for i in range(len(command)):
    print(command[i]);
    os.system(command[i]);
	
