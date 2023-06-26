using CSV, DataFrames, Plots, Plots.PlotMeasures, LaTeXStrings, StatsPlots
cd(@__DIR__)

raw = CSV.read("FA-sensitivity.csv",DataFrame);

struct MyStruct
    inventory
    salvage
    shortage
end


myvals = [];

for i in 1:nrow(raw)
    push!(
        myvals,
        MyStruct(
            raw."Inventory"[i],
            raw."Salvage"[i],
            raw."Shortage"[i])
        )
end



gb=groupby(raw, [:M, :penalty])

colnames = ["Inventory", "Salvage", "Shortage"]
ticklabel = [L"\nu = 0.05", L"\nu = 0.50", L"\nu = 5.00"]
penalty = [5,20,50];
ylb = Inf;
yub = 0;


for j in 1:2
    myplots = []
    for i in 1:3
        k = i + (j-1)*3
        if i == 1
            axval = true
            legnedval = true
            margval = 1mm
        else
            axval = :x
            legnedval = false
            margval = -15mm
        end
        if j == 2
            legnedval = false
        end
        vals = Matrix(gb[k][:,colnames]);
        ylb = min(ylb, minimum(vals))
        yub = max(yub, maximum(vals))
        p = groupedbar(
            vals,
            bar_position = :stack,
            bar_width=0.7,
            xticks=(1:3, ticklabel),
            label=reshape(colnames, (1,3)),
            xlabel=L"\textrm{penalty}~= %$(penalty[i])",
            legend=legnedval,
            left_margin=margval,
            #yticks = nothing,
            showaxis = axval,
            xtickfontsize = 11,
            ytickfontsize = 11,
            guidefont=font(15),
            legendfontsize = 15,
            gridalpha = 0.5,
            grid = :y,
            color = [:gray :blue :red]
            #minorticks = 5
        )
        push!(myplots,p)
    end
    layout = @layout [a b c];
    if j == 1
        #legnedval = false
        mytitle = L"\textrm{Emergency~procurement~allowed}"
    else
        #legnedval = true
        mytitle = L"\textrm{No~emergency~procurement~allowed}"
    end 

    p = plot(
        myplots[1], myplots[2], myplots[3], layout = layout,
        size = (1000,500),
        bottom_margin=5mm,
        top_margin = 0mm,
        link = :y,
        title = ["" mytitle ""],
        titlefontsize=20,
        )
    savefig(p, "temp$j.png")
end

