import os
import time
resultpath = "output/benchmark/"

command = []

instanceAttr3 = [" 0.0 "];
option = [" 0 ", " 1 "];

for k4 in range(len(option)):
    for k1 in range(len(instanceAttr3)):
        counter = 49+k1+100*k4;
        command = command + ["julia MAIN-case-study.jl " + str(counter) + instanceAttr3[k1] + " 1000 " + option[k4]];	
	

for i in range(len(command)):
    print(command[i]);
    os.system(command[i]);
	
