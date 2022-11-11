import os
import time
resultpath = "output/benchmark/"

command = []


command = command + ["julia MAIN.jl 2 3 10 0.1 1000 65 100 "];
command = command + ["julia MAIN.jl 4 3 10 0.2 1000 65 100 "];
command = command + ["julia MAIN.jl 6 3 10 0.3 1000 65 100 "];
command = command + ["julia MAIN.jl 8 3 10 0.4 1000 65 100 "];
command = command + ["julia MAIN.jl 10 3 10 0.5 1000 65 100 "];


for i in range(len(command)):
    print(command[i]);
    os.system(command[i]);
	
