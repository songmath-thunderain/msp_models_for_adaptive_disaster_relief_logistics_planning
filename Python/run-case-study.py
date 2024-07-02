import os
import time

command = []

instance_option = [" 2 "];
inv_level = [" 0 ", " 0.5 ", " 1 ", " 2 "]
fc_level = [" 0.25 ", " 0.5 ", " 1 ", " 2 ", " 4 "]
st_level = [" 3 ", " 4 "]

for k1 in instance_option:
    for k2 in inv_level:
            command = command + ["python main.py -p solveParams.yaml -d 0 -a 0 " + "-k 1 -o 1000 -ni 3 -nj 10 -t 1000000 " + "-s 1 -i " + k1 + "-w 1 -c 2 " + "-st 2 " + "-inv " + k2];

for k1 in instance_option:
    for k2 in fc_level:
        command = command + ["python main.py -p solveParams.yaml -d 0 -a 0 " + "-k 1 -o 1000 -ni 3 -nj 10 -t 1000000 " + "-s 1 -i " + k1 + "-w 1 -c 2 " + "-st 2 " + "-fc " + k2];


for k1 in instance_option:
    for k2 in st_level:
        command = command + ["python main.py -p solveParams.yaml -d 0 -a 0 " + "-k 1 -o 1000 -ni 3 -nj 10 -t 1000000 " + "-s 1 -i " + k1 + "-w 1 -c 2 " + "-st " + k2];



for i in range(len(command)):
    print(command[i]);
    os.system(command[i]);
	
