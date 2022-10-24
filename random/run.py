import os
import time
instancepath = "new_sample_instance/"
resultpath = "output/benchmark/"

command = []


instanceAttr1 = [" 3 ", " 6 ", " 9 "];
instanceAttr2 = [" 10 ", " 20 ", " 30 "];
instanceAttr3 = [" 0.005 ", " 0.05 ", " 0.5 ", " 5 "];

for k1 in range(len(instanceAttr3)):
    counter = 10*k1;
    for k2 in range(len(instanceAttr1)):
            counter = counter + 1;
            command = command + ["julia MAIN.jl " + str(counter) + instanceAttr1[k2] + instanceAttr2[k2] + instanceAttr3[k1] + " 1000 65 100 "];			

for i in range(len(command)):
    print(command[i]);
    os.system(command[i]);
	
