Random.seed!(2022);

max_iter = 100000;
stall = 200;
cutviol_maxiter = 100000;
nbhrs = 3;
time_limit = nbhrs*60^2;
ϵ = 1e-4;

########################################################################################
# Data input

# transition probability matrices:
P_intensity = Matrix(CSV.read("./case-study/mc_int_transition_prob.csv",DataFrame)) #intensity MC
Na = size(P_intensity)[1]; #intensity MC number of states

P_landfall = Matrix(CSV.read("./case-study/landfall_10.csv",DataFrame)) #landfall MC
Nc = size(P_landfall)[1]; #landfall MC number of states
T = copy(Nc); #define T_max
println("T = ", T);

# Now read in the track MC: note that each stage has their own track MC
trackMatrices = [];
trackStates = [];
for t = 1:(T-1)
	temp_name = "./case-study/mc_track_transition_prob_at_t"*string(t-1)*".csv";
	temp_MC = Matrix(CSV.read(temp_name,DataFrame));

	temp_name2 = "./case-study/mc_track_mean_error_at_t"*string(t)*".csv";
	temp_states = Matrix(CSV.read(temp_name2,DataFrame))[:,2];

	push!(trackMatrices, temp_MC);
	push!(trackStates, temp_states);

	if size(temp_MC)[2] != length(temp_states)
		print("Error in reading track MCs!");
		exit(0);
	end
end

########################################################################################
# Data processing into a joint MC

# First, count the total number of possible states
K = 1;
for t = 2:T
	global K += Na*length(trackStates[t-1]);
end

println("Total # of states K = ", K);

P_joint = zeros(K,K); #initialize the joint probability distribution MC
S = Array{Any,1}(undef,K); #list with elements [intensity,location-x,location-y], all are indices
absorbing_states = []; # list of absorbing states: Loc_y = T: landfall signifies the last stage (absorbing)
k1 = 1; #counter for the number of states
S[1] = [1,1,1];

for t = 2:T
	for l = 1:length(trackStates[t-1])
		for k = 1:Na
			global k1 += 1
			S[k1] = [k,l,t];
			if t == T
				push!(absorbing_states, k1);
			end
		end	
	end
end

println("absorbing_states = ", absorbing_states);

for k = 1:K
	if S[k][3] == T
		P_joint[k,k] = 1; # absorbing
	else
		for kk = 1:K
			if (S[k][2] <= size(trackMatrices[S[k][3]])[1]) && (S[kk][2] <= size(trackMatrices[S[k][3]])[2]) 
				P_joint[k,kk] = P_intensity[S[k][1],S[kk][1]]*trackMatrices[S[k][3]][S[k][2],S[kk][2]]*P_landfall[S[k][3],S[kk][3]];
			end
		end
	end
end

#normalize the probabilities
P_temp = deepcopy(P_joint);
for k=1:K, kk=1:K
    P_joint[k,kk] = P_temp[k,kk]/sum(P_temp[k,:]);
end

P_jointVec = deepcopy(P_joint);
# Get the smallest transition probability that is nonzero to give us a correct threshold to filter out impossible transitions
smallestTransProb = findmin(filter!(x -> x!= 0, vec(P_jointVec)))[1]*1.0/2; # JULIA WARNING: vec(P_joint) directly will change what P_joint looks like!!!!

println("smallestTransProb = ", smallestTransProb);


########################################################################################
# Data on the logistics network

d_JI = Matrix(CSV.read("./case-study/d_JI.csv",DataFrame));
d_II = Matrix(CSV.read("./case-study/d_II.csv",DataFrame));
d_KI = Matrix(CSV.read("./case-study/d_KI.csv",DataFrame));
d_KJ = Matrix(CSV.read("./case-study/d_KJ.csv",DataFrame));
d_SJ = Matrix(CSV.read("./case-study/d_SJ.csv",DataFrame));
x_cap = Matrix(CSV.read("./case-study/x_cap.csv",DataFrame));
max_D = Matrix(CSV.read("./case-study/demand.csv",DataFrame));

Nj = size(d_JI)[1]; # number of DPs
Ni = size(d_JI)[2]; # number of SPs
N0 = Ni + 1;

println("Ni = ", Ni, ", Nj = ", Nj);

#fuel cost
fuel = 0.0038;

#unit cost of transporting/rerouting items from MDC/SP i to/between SP i' 
cb = Array{Float64,3}(undef,N0,Ni,T);
for i=1:N0, ii=1:Ni, t=1:T
    if i < N0
        cb[i,ii,t] = fuel*d_II[i,ii]*(1+factor*(t-1))
    else
        cb[i,ii,t] = fuel*d_KI[ii]*(1+factor*(t-1))
    end 
end

#unit cost of transporting items from MDC/SP i to/between a demand point j
ca = Array{Float64,3}(undef,N0,Nj,T);
for i=1:N0, j=1:Nj, t=1:T
    if i < N0
        ca[i,j,t] = fuel*d_JI[j,i]*(1+factor*(t-1))
    else
        ca[i,j,t] = fuel*d_KJ[j]*(1+factor*(t-1))
    end 
end

base = 5; # base unit cost for logistic costs
h = Array{Float64,1}(undef,T); #unit cost for purchasing a relief item
ch = Array{Float64,2}(undef,Ni,T); #unit cost for holding an item at SP i
for t=1:T
    h[t] = base*(1+factor*(t-1));
	ch[:,t] = fill(0.05*base,Ni);
end

p = 50*base;
q = -0.2*base;

x_0 = zeros(Ni); #initial items at different SPs


########################################################################################
#Demand data 
x_max = 500;
max_intensity = 5;
SCEN = Array{Any,1}(undef,K); #SCEN is the list for each state k 

for k=1:K
    scen = zeros(Nj); #we will create a scenario list for each state k
    if S[k][3] == T
		#now lets calculate the demand from each DP 
		for j=1:Nj
			#how far did destination to DPj 
			dLandfall = d_SJ[S[k][2],j];
			if dLandfall >= x_max
				scen[j] = 0;
			else
				scen[j] = (1-dLandfall/x_max)*(S[k][1]/max_intensity)^2*max_D[j];
			end
		end
	end
    SCEN[k] = scen;
end

########################################################################################
########################################################################################
