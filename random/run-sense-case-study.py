import os
import time
resultpath = "output/benchmark/"

command = []

command = command + ["julia MAIN-case-study-sense.jl 101 0.05 1000 0 "];
command = command + ["julia MAIN-case-study-sense.jl 102 0.5 1000 0 "];
command = command + ["julia MAIN-case-study-sense.jl 103 5 1000 0 "];


command = command + ["julia MAIN-case-study-sense.jl 111 0.05 1000 1 "];
command = command + ["julia MAIN-case-study-sense.jl 112 0.5 1000 1 "];
command = command + ["julia MAIN-case-study-sense.jl 113 5 1000 1 "];

for i in range(len(command)):
    print(command[i]);
    os.system(command[i]);
	
