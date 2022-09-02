using StatsBase

function generate_out_of_sample(markov_chain::Array{Float64,2}, nb_samples::Int64, nb_stages::Int64, initial_state::Int64)
    out_of_sample = Array{Int64}(undef, nb_samples, nb_stages); 
    states = collect(1:size(markov_chain)[1])                   
    for s=1:nb_samples
        current_state = initial_state;
        for t=2:nb_stages
            out_of_sample[s,t-1] = current_state;
            current_state = sample(states, Weights(markov_chain[current_state,:]))
        end
        out_of_sample[s,nb_stages] = current_state;
    end
    return out_of_sample
end


# to test the function 
mc = [
    0.1 0.2  0.7;
    0.6 0.35 0.05;
    0.3 0.30 0.40];

out_of_sample = generate_out_of_sample(mc, 10, 5, 2);


# to save the file (in the current directory) 
using CSV, DataFrames
df = DataFrame(out_of_sample, :auto);
CSV.write("out_of_sample.csv", df)
pwd()