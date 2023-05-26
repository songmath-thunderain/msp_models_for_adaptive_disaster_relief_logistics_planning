import os
import time
resultpath = "output/benchmark/"

command = []

command = command + ["julia MAIN-sense.jl 1 6 10 0.05 1000 57 0 1 "];
command = command + ["julia MAIN-sense.jl 2 6 10 0.1 1000 57 0 1 "];
command = command + ["julia MAIN-sense.jl 3 6 10 0.15 1000 57 0 1 "];
command = command + ["julia MAIN-sense.jl 4 6 10 0.2 1000 57 0 1 "];
command = command + ["julia MAIN-sense.jl 5 6 10 0.25 1000 57 0 1 "];
command = command + ["julia MAIN-sense.jl 6 6 10 0.3 1000 57 0 1 "];
command = command + ["julia MAIN-sense.jl 7 6 10 0.35 1000 57 0 1 "];
command = command + ["julia MAIN-sense.jl 8 6 10 0.4 1000 57 0 1 "];
command = command + ["julia MAIN-sense.jl 9 6 10 0.45 1000 57 0 1 "];
command = command + ["julia MAIN-sense.jl 10 6 10 0.5 1000 57 0 1 "];

command = command + ["julia MAIN-sense.jl 11 6 10 0.05 1000 155 0 1 "];
command = command + ["julia MAIN-sense.jl 12 6 10 0.1 1000 155 0 1 "];
command = command + ["julia MAIN-sense.jl 13 6 10 0.15 1000 155 0 1 "];
command = command + ["julia MAIN-sense.jl 14 6 10 0.2 1000 155 0 1 "];
command = command + ["julia MAIN-sense.jl 15 6 10 0.25 1000 155 0 1 "];
command = command + ["julia MAIN-sense.jl 16 6 10 0.3 1000 155 0 1 "];
command = command + ["julia MAIN-sense.jl 17 6 10 0.35 1000 155 0 1 "];
command = command + ["julia MAIN-sense.jl 18 6 10 0.4 1000 155 0 1 "];
command = command + ["julia MAIN-sense.jl 19 6 10 0.45 1000 155 0 1 "];
command = command + ["julia MAIN-sense.jl 20 6 10 0.5 1000 155 0 1 "];

command = command + ["julia MAIN-sense.jl 21 6 10 0.05 1000 253 0 1 "];
command = command + ["julia MAIN-sense.jl 22 6 10 0.1 1000 253 0 1 "];
command = command + ["julia MAIN-sense.jl 23 6 10 0.15 1000 253 0 1 "];
command = command + ["julia MAIN-sense.jl 24 6 10 0.2 1000 253 0 1 "];
command = command + ["julia MAIN-sense.jl 25 6 10 0.25 1000 253 0 1 "];
command = command + ["julia MAIN-sense.jl 26 6 10 0.3 1000 253 0 1 "];
command = command + ["julia MAIN-sense.jl 27 6 10 0.35 1000 253 0 1 "];
command = command + ["julia MAIN-sense.jl 28 6 10 0.4 1000 253 0 1 "];
command = command + ["julia MAIN-sense.jl 29 6 10 0.45 1000 253 0 1 "];
command = command + ["julia MAIN-sense.jl 30 6 10 0.5 1000 253 0 1 "];

command = command + ["julia MAIN-sense.jl 31 6 10 0.05 1000 57 1 1 "];
command = command + ["julia MAIN-sense.jl 32 6 10 0.1 1000 57 1 1 "];
command = command + ["julia MAIN-sense.jl 33 6 10 0.15 1000 57 1 1 "];
command = command + ["julia MAIN-sense.jl 34 6 10 0.2 1000 57 1 1 "];
command = command + ["julia MAIN-sense.jl 35 6 10 0.25 1000 57 1 1 "];
command = command + ["julia MAIN-sense.jl 36 6 10 0.3 1000 57 1 1 "];
command = command + ["julia MAIN-sense.jl 37 6 10 0.35 1000 57 1 1 "];
command = command + ["julia MAIN-sense.jl 38 6 10 0.4 1000 57 1 1 "];
command = command + ["julia MAIN-sense.jl 39 6 10 0.45 1000 57 1 1 "];
command = command + ["julia MAIN-sense.jl 40 6 10 0.5 1000 57 1 1 "];

command = command + ["julia MAIN-sense.jl 41 6 10 0.05 1000 155 1 1 "];
command = command + ["julia MAIN-sense.jl 42 6 10 0.1 1000 155 1 1 "];
command = command + ["julia MAIN-sense.jl 43 6 10 0.15 1000 155 1 1 "];
command = command + ["julia MAIN-sense.jl 44 6 10 0.2 1000 155 1 1 "];
command = command + ["julia MAIN-sense.jl 45 6 10 0.25 1000 155 1 1 "];
command = command + ["julia MAIN-sense.jl 46 6 10 0.3 1000 155 1 1 "];
command = command + ["julia MAIN-sense.jl 47 6 10 0.35 1000 155 1 1 "];
command = command + ["julia MAIN-sense.jl 48 6 10 0.4 1000 155 1 1 "];
command = command + ["julia MAIN-sense.jl 49 6 10 0.45 1000 155 1 1 "];
command = command + ["julia MAIN-sense.jl 50 6 10 0.5 1000 155 1 1 "];

command = command + ["julia MAIN-sense.jl 51 6 10 0.05 1000 253 1 1 "];
command = command + ["julia MAIN-sense.jl 52 6 10 0.1 1000 253 1 1 "];
command = command + ["julia MAIN-sense.jl 53 6 10 0.15 1000 253 1 1 "];
command = command + ["julia MAIN-sense.jl 54 6 10 0.2 1000 253 1 1 "];
command = command + ["julia MAIN-sense.jl 55 6 10 0.25 1000 253 1 1 "];
command = command + ["julia MAIN-sense.jl 56 6 10 0.3 1000 253 1 1 "];
command = command + ["julia MAIN-sense.jl 57 6 10 0.35 1000 253 1 1 "];
command = command + ["julia MAIN-sense.jl 58 6 10 0.4 1000 253 1 1 "];
command = command + ["julia MAIN-sense.jl 59 6 10 0.45 1000 253 1 1 "];
command = command + ["julia MAIN-sense.jl 60 6 10 0.5 1000 253 1 1 "];

for i in range(len(command)):
    print(command[i]);
    os.system(command[i]);
	
