import os
import time

command = []

instance_option = [" 1 "];

st_levels = [" 1 ", " 3 ", " 4 "];

fc_levels = [" 0.25 ", " 0.5 ", " 2 ", " 4 "]

for k1 in instance_option:
    for k3 in st_levels:
        command = command + ["python main.py -p solveParams.yaml -d 0 -a 0 " + "-k 1 -o 1000 -ni 3 -nj 10 -t 2 " + "-s 5 -i " + k1 + "-w 0 -c 1 " + "-st " + k3 + "-trans 10 "];

    for k3 in fc_levels:
        command = command + ["python main.py -p solveParams.yaml -d 0 -a 0 " + "-k 1 -o 1000 -ni 3 -nj 10 -t 2 " + "-s 5 -i " + k1 + "-w 0 -c 1 " + "-st 2 " + "-fc " + k3 + "-trans 10 "];

for i in range(len(command)):
    print(command[i]);
    os.system(command[i]);
	
