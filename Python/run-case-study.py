import os
import time

command = []

instance_option = [" 2 "];
cost_structure = [" 1 "];
tau = [" 5 ", " 10 "]
st = [" 3 ", " 4 "]

for k1 in instance_option:
    for k2 in tau:
        for k3 in st:
            command = command + ["python main.py -p solveParams.yaml -d 0 -a 0 " + "-k 1 -o 1000 -ni 3 -nj 10 -t " + k2 + "-s 5 -i " + k1 + "-w 1 -c 1 " + "-st " + k3];

for i in range(len(command)):
    print(command[i]);
    os.system(command[i]);
	
