import os
import time
resultpath = "output/benchmark/"

command = []
"""
command = command + ["julia MAIN-case-study-sense.jl 101 0.05 1000 0 5 "];
command = command + ["julia MAIN-case-study-sense.jl 102 0.5 1000 0 5 "];
command = command + ["julia MAIN-case-study-sense.jl 103 5 1000 0 5 "];


command = command + ["julia MAIN-case-study-sense.jl 111 0.05 1000 1 5 "];
command = command + ["julia MAIN-case-study-sense.jl 112 0.5 1000 1 5 "];
command = command + ["julia MAIN-case-study-sense.jl 113 5 1000 1 5 "];

command = command + ["julia MAIN-case-study-sense.jl 121 0.05 1000 0 20 "];
command = command + ["julia MAIN-case-study-sense.jl 122 0.5 1000 0 20 "];
command = command + ["julia MAIN-case-study-sense.jl 123 5 1000 0 20 "];


command = command + ["julia MAIN-case-study-sense.jl 131 0.05 1000 1 20 "];
command = command + ["julia MAIN-case-study-sense.jl 132 0.5 1000 1 20 "];
command = command + ["julia MAIN-case-study-sense.jl 133 5 1000 1 20 "];
"""
command = command + ["julia MAIN-case-study-sense.jl 161 0.05 1000 0 100 "];
command = command + ["julia MAIN-case-study-sense.jl 162 0.5 1000 0 100 "];
command = command + ["julia MAIN-case-study-sense.jl 163 5 1000 0 100 "];


command = command + ["julia MAIN-case-study-sense.jl 171 0.05 1000 1 100 "];
command = command + ["julia MAIN-case-study-sense.jl 172 0.5 1000 1 100 "];
command = command + ["julia MAIN-case-study-sense.jl 173 5 1000 1 100 "];

for i in range(len(command)):
    print(command[i]);
    os.system(command[i]);
	
