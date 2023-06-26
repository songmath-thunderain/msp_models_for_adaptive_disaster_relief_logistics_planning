using CSV, DataFrames, Plots, Plots.PlotMeasures, LaTeXStrings, XLSX#, GR

cd(@__DIR__)
xf = XLSX.readxlsx("procurement-percentage-case-study.xlsx")
sh = xf["all"] # get a reference to a Worksheet
dt = sh["A1:I3"];
rhos = [1,2,3,4,5,6,7,8,9];
vals = convert.(Float64,dt);
vals = transpose(vals);
policies = [L"\nu = 0.05", L"\nu = 0.5", L"\nu = 5"];

#p = gr()
p = plot(
    rhos,
    vals,
    xlabel=L"t",
    ylabel=L"\textrm{Proc}~(\%)",
    yformatter=y->(string(floor(Int,y*100))*"%"),
    label=[policies[1] policies[2] policies[3]],
    #lw=7,
    lw=3,
    linecolor = [:black :blue :red],
    #linestyle=[:solid :dash :dashdotdot],
    #linecolor=[:black :black :black],
    #marker = [:square :circle :diamond],
    #markercolor = [:black :black :black],
    #markersize = [5 5 5],
    #color = [:black :gray :blue :red],
    #leg=(0.09,0.955),
    leg=false,
    #windowsize=(1000,550),
    #windowsize=(500,250),
    xtickfontsize = 11,
    ytickfontsize = 11,
    guidefont=font(25),
    #yguidefontrotation=-90,
    yguideposition = :left,
    yguidefontvalign=:center,
    #yticks=7,
    #ylims=collect(1:maximum(avg_data)/10:maximum(avg_data)),
    #ylims=(0,maximum(avg_data)),
    #bottom_margin=2mm,
    left_margin=9mm,
    #frame=:box,
    legendfontsize = 13,
    widen = true,
    fmt = :png,
    minorgrid = true,
    minorticks = 2,
        )
savefig(p, "procurement-percentage-case-study.png")
