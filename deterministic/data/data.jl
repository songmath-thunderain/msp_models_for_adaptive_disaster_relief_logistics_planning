#fix the random seed
Random.seed!(2022);

#Ni: number of supply points (without MDC).
#Nj: number of demand points.
#N0: number of supply points including MDC.
#x_cap: the capacity of each supply point.
#cb: unit cost of transporting/rerouting relief items from the MDC/SP i to/between different SP i' 
#ch: unit cost for holding an item at SP i.
#h: unit cost for purchasing a relief item.
#-----------------------------------------------------------------------------------------
#ca: unit cost of transporting relief items from the MDC/SP i to/between a demand point j.
#p: penalty cost for failing to serve the demand.
#q: salvage cost for each unit of overstock.
#-----------------------------------------------------------------------------------------
#Na: number of states in the intensity MC
#Nb: number of states in the locations MC
#K: number of states in the joint distrubutions between intesities and locations
#P_intesity: transition probability matrix of intensity MC
#P_location: transition probability matrix of location MC
#P_joint: transition probability matrix of the joint distrubution beween intensity and location MC
#M: is the number of samples in the second layer [REVISION: removed]
#-----------------------------------------------------------------------------------------
#right now the storage capacities are uniformly distributed between (0.01,1.18) max_demand
#if we were to add a restriction on the flow, it should be unifromly distributed between (0.05,4.6)

#M = 10; #number of samples in the second layer [REVISION]
max_iter = 100000;
stall = 200;
cutviol_maxiter = 100000;
nbhrs = 3;
time_limit = nbhrs*60^2;
ϵ = 1e-4;

########################################################################################
########################################################################################


N0=Ni+1; #number of supply points + the MDC

#the set of possible locations
L = [(0,100),(100,200),(200,300),(300,400),(400,500),(500,600),(600,700)]; #discussion on the last state

# probability distributions:
P_intesity = Matrix(CSV.read("./data/intensity.csv",DataFrame)) #intensity MC
P_location = Matrix(CSV.read("./data/location.csv",DataFrame)) #location MC
#P_landfall = Matrix(CSV.read(string("./data/landfall_",Tmax,".csv"),DataFrame)) #landfall MC

Na = size(P_intesity)[1]; #intensity MC number of states
Nb = size(P_location)[1]; #location MC number of states
#Nc = size(P_landfall)[1]; #landfall MC number of states

#exp_Ts = expected_number_of_stages(P_landfall)
T =  Tmax; #Note: Tmax is user-defined [REVISION]
K = Na*Nb;

P_joint = zeros(K,K); #initialize the joint probability distribution between the intensity and location MC
S = Array{Any,1}(undef,K); #a list to store what each state correspond to, with its elements being [intensity,location]

#initialState = (initialNa-1)*Na + initialNb;
#initialState = OS_paths[1,1];
#create the joint transition probability distribution matrix 
k1 = 0;
for k=1:Na, l=1:Nb
    global k1 +=1
    k2 = 0;
    for n=1:Na, m=1:Nb
        k2 +=1
        P_joint[k1,k2] = P_intesity[k,n]*P_location[l,m];
    end
    S[k1] = [k,l]; 
end

# Create the transition probability from stage 1 to stage T (applying the C-K equation)
P_terminals = Matrix{Float64}[]
for t = 1:T  
    P_terminal = P_joint^(T-t);
    P_terminal_c = deepcopy(P_terminal);
    #normalize the probabilities
    for k=1:K, kk=1:K
       P_terminal[k,kk] = P_terminal_c[k,kk]/sum(P_terminal_c[k,:])
    end
    push!(P_terminals, P_terminal);
end    


#normalize the probabilities
P_temp = deepcopy(P_joint);
for k=1:K, kk=1:K
    P_joint[k,kk] = P_temp[k,kk]/sum(P_temp[k,:])
end
    
########################################################################################

#locations of the different facilities:
MDC = [350,450];

x_low = 0;
x_up = 700;

y_low = 0;
y_up = 100;

nodes = CSV.read("./data/nodes.csv",DataFrame);

SP = []; #initialize a list for the coordinates of the different supply points
DP = []; #initialize a list for the coordinates of the different demand points

for i=1:Ni
    push!(SP,[nodes[i,1],nodes[i,2]+100]);
end

for j=1:Nj
    push!(DP,[nodes[j,3],nodes[j,4]]);
end

#fuel cost
fuel = 0.0038;


#unit cost of transporting/rerouting relief items from the MDC/SP i to/between different SP i' 
cb = Array{Float64,3}(undef,N0,Ni,T);
for i=1:N0, ii=1:Ni, t=1:T
    if i < N0
        cb[i,ii,t] = fuel*norm(SP[i]-SP[ii],2)*(1+factor*(t-1))
    else
        cb[i,ii,t] = fuel*norm(MDC-SP[ii],2)*(1+factor*(t-1))
    end 
end

#unit cost of transporting relief items from the MDC/SP i to/between a demand point j
ca = Array{Float64,3}(undef,N0,Nj,T);
for i=1:N0, j=1:Nj, t=1:T
    if i < N0
        ca[i,j,t] = fuel*norm(SP[i]-DP[j],2)*(1+factor*(t-1))
    else
        ca[i,j,t] = fuel*norm(MDC-DP[j],2)*(1+factor*(t-1))
    end 
end


h = Array{Float64,1}(undef,T);
ch = Array{Float64,2}(undef,Ni,T);

# base unit cost for purchasing a relief item
base = 5;

#p: penalty cost for failing to serve the demand
#p = 80*base; 
p = 50*base;

#q: salvage cost for each unit of overstock.
#q = -0.05*base;
q = -0.2*base;

for t=1:T
    #h: unit cost for purchasing a relief item
    h[t] = base*(1+factor*(t-1));
    
    #ch: unit cost for holding an item at SP i (2% of the purchasing cost)
    ch[:,t] = fill(0.05*base,Ni);
	#ch[:,t] = fill(0.2*base,Ni);
end

#the maximum demand that can happen
D_max = 400;

# [REVISION]: to be consistent with random 
x_cap = Matrix(nodes)[1:Ni,5]*(Nj/Ni); #capacity of each SP
x_0 = zeros(Ni); #initial items at different SPs
f_cap = fill(Inf,N0,Ni);

########################################################################################
#Demand data 

#Note: we get rid of layer2 in the revised code! [REVISION]
#=
layer2 = [];

for l in L
    cover = l;
    list = [];

    #left = [-cover[2],-cover[1]]
    right = [cover[1],cover[2]]

    #parts_left = (left[2]-left[1])/M;
    parts_right = (right[2]-right[1])/(2*M);

    count_right = right[1];
    count_2 = 0
    while true    
        count_2+=1
        count_right +=parts_right
        if count_right >= right[2]
            break
        elseif isodd(count_2)
            push!(list,count_right)
        end
    end
    push!(layer2,list);
end
=#

#SCEN is the list for each state k 
SCEN = Array{Any,1}(undef,K);

#the largest possible radius of a hurricane
c_max = 300;

#=
# Note: we get rid of layer2 in the revised code! [REVISION]
for k=1:K
    #initialize a scenario matrix (scen) for each DP j and each possible point m in layer2;
    #we will create a scenario list for each state k
    scen = zeros(Nj,M);
    
    a = S[k][1]; #what is the observed intensity state
    l = S[k][2]; #what is the observed location state
    
    #what are the coordinates of where the hurricane made landfall [xx,yy]
    #note that yy=0 since the huricane makes landfall in the cost by assumption
    #since we observed state location l, here are the M possible locations
    for m=1:M
        predicted = layer2[l][m]; #this is the predicted x_coordinates for the landfall location
        xx_coord = max(x_low,min(predicted,x_up)); #we already now x in can't be smaller than x_low and bigger than x_up; 
        landfall = [xx_coord,0]
        #now lets calculate the demand from each DP to the m location 
        for j=1:Nj
            #how far did destination to DPj j
            c_j = norm(landfall-DP[j],2);
            if c_j <= c_max
                scen[j,m] = D_max*(1-(c_j/c_max))*(a-1)^2/((Na-1)^2)
            else
                scen[j,m] = 0;
            end
        end
    end
    SCEN[k] = scen;
end
=#

for k=1:K
    #initialize a scenario matrix (scen) for each DP j and each possible landfall location;
    #we will create a scenario list for each state k
    scen = zeros(Nj);
    
    a = S[k][1]; #what is the observed intensity state
    l = S[k][2]; #what is the observed location state
    
    #what are the coordinates of where the hurricane made landfall [xx,yy]
    #note that yy=0 since the huricane makes landfall in the cost by assumption
	predicted = 0.5*(L[l][1]+L[l][2]); #this is the predicted x_coordinates for the landfall location: for simplicity, just use the center of the interval
	xx_coord = max(x_low,min(predicted,x_up)); #we already now x in can't be smaller than x_low and bigger than x_up; 
	landfall = [xx_coord,0]
	#now lets calculate the demand from each DP to the location 
	for j=1:Nj
		#how far did destination to DPj j
		c_j = norm(landfall-DP[j],2);
		if c_j <= c_max
			scen[j] = D_max*(1-(c_j/c_max))*(a-1)^2/((Na-1)^2)
		else
			scen[j] = 0;
		end
	end
    SCEN[k] = scen;
end


########################################################################################
########################################################################################
