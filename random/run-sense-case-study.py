import os
import time
resultpath = "output/benchmark/"

command = []

command = command + ["julia MAIN-case-study-sense.jl 101 0.05 1000 0 2 "];
command = command + ["julia MAIN-case-study-sense.jl 102 0.5 1000 0 2 "];
command = command + ["julia MAIN-case-study-sense.jl 103 5 1000 0 2 "];


command = command + ["julia MAIN-case-study-sense.jl 111 0.05 1000 1 2 "];
command = command + ["julia MAIN-case-study-sense.jl 112 0.5 1000 1 2 "];
command = command + ["julia MAIN-case-study-sense.jl 113 5 1000 1 2 "];

command = command + ["julia MAIN-case-study-sense.jl 121 0.05 1000 0 10 "];
command = command + ["julia MAIN-case-study-sense.jl 122 0.5 1000 0 10 "];
command = command + ["julia MAIN-case-study-sense.jl 123 5 1000 0 10 "];


command = command + ["julia MAIN-case-study-sense.jl 131 0.05 1000 1 10 "];
command = command + ["julia MAIN-case-study-sense.jl 132 0.5 1000 1 10 "];
command = command + ["julia MAIN-case-study-sense.jl 133 5 1000 1 10 "];

command = command + ["julia MAIN-case-study-sense.jl 141 0.05 1000 0 50 "];
command = command + ["julia MAIN-case-study-sense.jl 142 0.5 1000 0 50 "];
command = command + ["julia MAIN-case-study-sense.jl 143 5 1000 0 50 "];


command = command + ["julia MAIN-case-study-sense.jl 151 0.05 1000 1 50 "];
command = command + ["julia MAIN-case-study-sense.jl 152 0.5 1000 1 50 "];
command = command + ["julia MAIN-case-study-sense.jl 153 5 1000 1 50 "];

for i in range(len(command)):
    print(command[i]);
    os.system(command[i]);
	
