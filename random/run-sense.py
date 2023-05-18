import os
import time
resultpath = "output/benchmark/"

command = []

command = command + ["julia MAIN.jl 1 9 30 0.05 1000 57 0 1 "];
command = command + ["julia MAIN.jl 2 9 30 0.1 1000 57 0 1 "];
command = command + ["julia MAIN.jl 3 9 30 0.15 1000 57 0 1 "];
command = command + ["julia MAIN.jl 4 9 30 0.2 1000 57 0 1 "];
command = command + ["julia MAIN.jl 5 9 30 0.25 1000 57 0 1 "];
command = command + ["julia MAIN.jl 6 9 30 0.3 1000 57 0 1 "];
command = command + ["julia MAIN.jl 7 9 30 0.35 1000 57 0 1 "];
command = command + ["julia MAIN.jl 8 9 30 0.4 1000 57 0 1 "];
command = command + ["julia MAIN.jl 9 9 30 0.45 1000 57 0 1 "];
command = command + ["julia MAIN.jl 10 9 30 0.5 1000 57 0 1 "];

command = command + ["julia MAIN.jl 11 9 30 0.05 1000 155 0 1 "];
command = command + ["julia MAIN.jl 12 9 30 0.1 1000 155 0 1 "];
command = command + ["julia MAIN.jl 13 9 30 0.15 1000 155 0 1 "];
command = command + ["julia MAIN.jl 14 9 30 0.2 1000 155 0 1 "];
command = command + ["julia MAIN.jl 15 9 30 0.25 1000 155 0 1 "];
command = command + ["julia MAIN.jl 16 9 30 0.3 1000 155 0 1 "];
command = command + ["julia MAIN.jl 17 9 30 0.35 1000 155 0 1 "];
command = command + ["julia MAIN.jl 18 9 30 0.4 1000 155 0 1 "];
command = command + ["julia MAIN.jl 19 9 30 0.45 1000 155 0 1 "];
command = command + ["julia MAIN.jl 20 9 30 0.5 1000 155 0 1 "];

command = command + ["julia MAIN.jl 21 9 30 0.05 1000 253 0 1 "];
command = command + ["julia MAIN.jl 22 9 30 0.1 1000 253 0 1 "];
command = command + ["julia MAIN.jl 23 9 30 0.15 1000 253 0 1 "];
command = command + ["julia MAIN.jl 24 9 30 0.2 1000 253 0 1 "];
command = command + ["julia MAIN.jl 25 9 30 0.25 1000 253 0 1 "];
command = command + ["julia MAIN.jl 26 9 30 0.3 1000 253 0 1 "];
command = command + ["julia MAIN.jl 27 9 30 0.35 1000 253 0 1 "];
command = command + ["julia MAIN.jl 28 9 30 0.4 1000 253 0 1 "];
command = command + ["julia MAIN.jl 29 9 30 0.45 1000 253 0 1 "];
command = command + ["julia MAIN.jl 30 9 30 0.5 1000 253 0 1 "];

command = command + ["julia MAIN.jl 31 9 30 0.05 1000 57 1 1 "];
command = command + ["julia MAIN.jl 32 9 30 0.1 1000 57 1 1 "];
command = command + ["julia MAIN.jl 33 9 30 0.15 1000 57 1 1 "];
command = command + ["julia MAIN.jl 34 9 30 0.2 1000 57 1 1 "];
command = command + ["julia MAIN.jl 35 9 30 0.25 1000 57 1 1 "];
command = command + ["julia MAIN.jl 36 9 30 0.3 1000 57 1 1 "];
command = command + ["julia MAIN.jl 37 9 30 0.35 1000 57 1 1 "];
command = command + ["julia MAIN.jl 38 9 30 0.4 1000 57 1 1 "];
command = command + ["julia MAIN.jl 39 9 30 0.45 1000 57 1 1 "];
command = command + ["julia MAIN.jl 40 9 30 0.5 1000 57 1 1 "];

command = command + ["julia MAIN.jl 41 9 30 0.05 1000 155 1 1 "];
command = command + ["julia MAIN.jl 42 9 30 0.1 1000 155 1 1 "];
command = command + ["julia MAIN.jl 43 9 30 0.15 1000 155 1 1 "];
command = command + ["julia MAIN.jl 44 9 30 0.2 1000 155 1 1 "];
command = command + ["julia MAIN.jl 45 9 30 0.25 1000 155 1 1 "];
command = command + ["julia MAIN.jl 46 9 30 0.3 1000 155 1 1 "];
command = command + ["julia MAIN.jl 47 9 30 0.35 1000 155 1 1 "];
command = command + ["julia MAIN.jl 48 9 30 0.4 1000 155 1 1 "];
command = command + ["julia MAIN.jl 49 9 30 0.45 1000 155 1 1 "];
command = command + ["julia MAIN.jl 50 9 30 0.5 1000 155 1 1 "];

command = command + ["julia MAIN.jl 51 9 30 0.05 1000 253 1 1 "];
command = command + ["julia MAIN.jl 52 9 30 0.1 1000 253 1 1 "];
command = command + ["julia MAIN.jl 53 9 30 0.15 1000 253 1 1 "];
command = command + ["julia MAIN.jl 54 9 30 0.2 1000 253 1 1 "];
command = command + ["julia MAIN.jl 55 9 30 0.25 1000 253 1 1 "];
command = command + ["julia MAIN.jl 56 9 30 0.3 1000 253 1 1 "];
command = command + ["julia MAIN.jl 57 9 30 0.35 1000 253 1 1 "];
command = command + ["julia MAIN.jl 58 9 30 0.4 1000 253 1 1 "];
command = command + ["julia MAIN.jl 59 9 30 0.45 1000 253 1 1 "];
command = command + ["julia MAIN.jl 60 9 30 0.5 1000 253 1 1 "];

for i in range(len(command)):
    print(command[i]);
    os.system(command[i]);
	
