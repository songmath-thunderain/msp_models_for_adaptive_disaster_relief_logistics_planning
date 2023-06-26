using CSV, DataFrames, Plots, Plots.PlotMeasures, StatsPlots, LaTeXStrings, Statistics
cd(@__DIR__)

raw = CSV.read("data_to_convert_to_table3.csv",DataFrame);
data = Matrix(raw);

df = DataFrame("FA-MSP" => data[:,4] , "Wait-and-see" => data[:,8], "S-2SSP" => data[:,11], "RH-2SSP" => data[:,15]);
matrix = Matrix(df);

rho = ["0.05", "0.5", "5"];
avg_data = Array{Float64, 2}(undef, 3, 4);
for r=1:length(rho)
    from = (r-1)*9+1
    to = (r-1)*9+9
    avg_data[r,:] = [mean(matrix[from:to,k]) for k=1:4] 
end


policies = ["FA-MSP","Decision-tree","RH-2SSP","S-2SSP"];

p = groupedbar(
    rho,
    avg_data,
    xlabel=L"\nu",
    ylabel=L"\textrm{Gap}~(\%)~\textrm{on~average}",
    yformatter=y->(string(floor(Int,y))*"%"),
    label=[policies[1] policies[2] policies[3] policies[4]],
    #bar_position = :dodge,
    #bar_width=0.7,
    color = [:black :gray :blue :red],
    #leg=(0.09,0.955),
    leg=:topright,
    #windowsize=(1000,550),
    xtickfontsize = 10,
    ytickfontsize = 10,
    guidefont=font(20),
    #yticks=7,
    #ylims=collect(1:maximum(avg_data)/10:maximum(avg_data)),
    ylims=(0,maximum(avg_data)),
    #bottom_margin=10mm,
    #left_margin=5mm,
    #frame=:box,
    legendfontsize = 10,
    #widen = true,
    fmt = :png
    )
savefig(p, "all_avg.png")

