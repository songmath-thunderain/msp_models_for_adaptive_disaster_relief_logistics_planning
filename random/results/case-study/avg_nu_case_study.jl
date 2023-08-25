using CSV, DataFrames, Plots, Plots.PlotMeasures, StatsPlots, LaTeXStrings, Statistics

rho = ["0.05", "0.5", "5"];
avg_data = Array{Float64, 2}(undef, 3, 3);
avg_data[1,:] = [3,135,165]
avg_data[2,:] = [73,117,154]
avg_data[3,:] = [82,82,90]

policies = ["FA-MSP","RH-2SSP","S-2SSP"];

p = groupedbar(
    rho,
    avg_data,
    xlabel=L"\nu",
    ylabel=L"\textrm{Gap}~(\%)~\textrm{on~average}",
    label=[policies[1] policies[2] policies[3]],
    #bar_position = :dodge,
    #bar_width=0.7,
    color = [:black :blue :red],
    #leg=(0.09,0.955),
    leg=:topright,
    #windowsize=(1000,550),
    xtickfontsize = 10,
    ytickfontsize = 10,
    guidefont=font(20),
    yticks=7,
    #ylims=collect(1:maximum(avg_data)/10:maximum(avg_data)),
    ylims=(0,maximum(avg_data)),
    #bottom_margin=10mm,
    #left_margin=5mm,
    #frame=:box,
    legendfontsize = 10,
    #widen = true,
    fmt = :png
    )
savefig(p, "all_case_study_option1.png")
