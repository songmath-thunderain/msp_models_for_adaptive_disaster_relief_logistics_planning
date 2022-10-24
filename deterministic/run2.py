import os
import time
instancepath = "new_sample_instance/"
resultpath = "output/benchmark/"

command = []

command = command + ["julia MAINCV.jl 15 6  20  0.05  1000 9 7 "];
command = command + ["julia MAINCV.jl 21 3  10  0.5  1000 9 7 "];
command = command + ["julia MAINCV.jl 35 6  20  5  1000 9 7 "];

for i in range(len(command)):
    print(command[i]);
    os.system(command[i]);
	

