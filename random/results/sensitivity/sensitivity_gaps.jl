using CSV, DataFrames, Plots, Plots.PlotMeasures, LaTeXStrings, XLSX#, GR

cd(@__DIR__)
myplots= []
for m in [1,3]

    xf = XLSX.readxlsx("all$(m).xlsx")
    sh = xf["all"] # get a reference to a Worksheet
    dt = sh["F1:J11"];
    rhos = convert.(Float64,dt[2:size(dt)[1],1]);
    vals = convert.(Float64,dt[2:size(dt)[1],2:end]);
    policies = dt[1,2:end];
    policies = ["FA-MSP", "Decision-tree", "RH-2SSP", "S-2SSP"];

    #labelflag=(0.835,0.9825)
    labelflag=false
    if m == 1
        #labelflag=:right
        labelflag=(0.15,0.55)
    end
    #p = gr()
    p = plot(
        rhos,
        vals,
        xlabel=L"\nu",
        ylabel= L"\textrm{Gap}~(\%)",
        yformatter=y->(string(floor(Int,y*100))*"%"),
        label=[policies[1] policies[2] policies[3] policies[4]],
        legend=labelflag,
        lw=3,
        linecolor = [:black :gray :blue :red],
        #linestyle=[:solid :dash :dashdotdot :dot],
        #linecolor=[:black :black :black],
        #marker = [:square :circle :diamond],
        #markercolor = [:black :black :black],
        #markersize = [5 5 5],
        #color = [:black :gray :blue :red],
        #leg=(0.835,0.9825),
        #leg=:topright,
        windowsize=(500,300),
        xtickfontsize = 10,
        ytickfontsize = 10,
        guidefont=font(20),
        #yticks=7,
        #ylims=collect(1:maximum(avg_data)/10:maximum(avg_data)),
        #ylims=(0,maximum(avg_data)),
        #bottom_margin=10mm,
        left_margin=5mm,
        #right_margin=5mm,
        #frame=:box,
        legendfontsize = 13,
        #widen = true,
        fmt = :png
            )
    push!(myplots,p)
    savefig(p, "gaps$(m).png")
end
myplots[1]