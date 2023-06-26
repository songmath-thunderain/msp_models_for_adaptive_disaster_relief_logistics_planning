using CSV, DataFrames, Plots, Plots.PlotMeasures, LaTeXStrings, XLSX#, GR

cd(@__DIR__)
xf = XLSX.readxlsx("procurement-avgAmount-case-study.xlsx")
sh = xf["all"] # get a reference to a Worksheet
dt = sh["A1:I3"];
rhos = [1,2,3,4,5,6,7,8,9];
vals = convert.(Float64,dt);
vals = transpose(vals);
policies = [L"~\nu = 0.05", L"~\nu = 0.5", L"~\nu = 5"];

#p = gr()
p = plot(
    rhos,
    vals,
    xlabel=L"t",
    ylabel=L"\bar{f}_t",
    #yformatter=y->(string(floor(Int,y*100))*"%"),
    label=[policies[1] policies[2] policies[3]],
    lw=3,
    linecolor = [:black :blue :red],
    #linestyle=[:solid :dot :dashdotdot],
    leg=:topright,
    #windowsize=(500,250),
    #xtickfontsize = 10,
    #ytickfontsize = 10,
    guidefont=font(20),
    yguidefontrotation=-90,
    yguideposition = :left,
    yguidefontvalign=:center,
    #yticks=7,
    #ylims=collect(1:maximum(avg_data)/10:maximum(avg_data)),
    #ylims=(0,maximum(avg_data)),
    #bottom_margin=2mm,
    left_margin=8mm,
    #frame=:box,
    legendfontsize = 12,
    widen = true,
    fmt = :png,
    minorgrid = true,
    minorticks = 2
        )
savefig(p, "procurement-avgAmount-case-study.png")
p

p = plot(
    rhos,
    vals,
    xlabel=L"t",
    ylabel=L"\bar{f}_t",
    #yformatter=y->(string(floor(Int,y*100))*"%"),
    label=[policies[1] policies[2] policies[3]],
    lw=3,
    linecolor = [:black :blue :red],
    #linestyle=[:solid :dot :dashdotdot],
    leg=:top,
    #windowsize=(500,250),
    xtickfontsize = 11,
    ytickfontsize = 11,
    guidefont=font(25),
    yguidefontrotation=-90,
    yguideposition = :left,
    yguidefontvalign=:center,
    #yticks=7,
    #ylims=collect(1:maximum(avg_data)/10:maximum(avg_data)),
    #ylims=(0,maximum(avg_data)),
    #bottom_margin=2mm,
    left_margin=9mm,
    #frame=:box,
    legendfontsize = 12,
    widen = true,
    fmt = :png,
    minorgrid = true,
    minorticks = 2,
    legend_column = -1
        )
#lens!([1, 6], [0.9, 1.1], inset = (1, bbox(0.5, 0.0, 0.4, 0.4)))
p

lens!([1, 9], [0.0, 3000], inset = (1, bbox(0.3, 0.25, 0.7, 0.4)))
savefig(p, "procurement-avgAmount-case-study.png")