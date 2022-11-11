import os
import time
resultpath = "output/benchmark/"

command = []


instanceAttr1 = [" 177 ", " 289 "];
instanceAttr3 = [" 0.05 ", " 0.1 ", " 0.15 ", "0.2", " 0.25 ", " 0.3 ", " 0.35 ", "0.4", " 0.45 ", " 0.5 "];

for k1 in range(len(instanceAttr1)):
    counter = 12*(k1+1);
    for k3 in range(len(instanceAttr3)):
            counter = counter + 1;
            command = command + ["julia MAIN.jl " + str(counter) + " 3 " + " 10 " + instanceAttr3[k3] + " 1000 " + instanceAttr1[k1] + " 100 "];			

for i in range(len(command)):
    print(command[i]);
    os.system(command[i]);
	
