import os
import time
resultpath = "output/benchmark/"

command = []

instanceAttr1 = [" 3 ", " 6 ", " 9 "];
instanceAttr2 = [" 10 ", " 20 ", " 30 "];
instanceAttr3 = [" 0.05 ", " 0.5 ", " 5 "];
option = [" 0 ", " 1 "];

for k4 in range(len(option)):
    for k1 in range(len(instanceAttr3)):
        counter = 10*k1+100*k4;
        for k2 in range(len(instanceAttr1)):
            for k3 in range(len(instanceAttr2)):
                counter = counter + 1;
                command = command + ["julia MAIN.jl " + str(counter) + instanceAttr1[k2] + instanceAttr2[k3] + instanceAttr3[k1] + " 1000 57 " + option[k4]];	
	

for i in range(len(command)):
    print(command[i]);
    os.system(command[i]);
	
