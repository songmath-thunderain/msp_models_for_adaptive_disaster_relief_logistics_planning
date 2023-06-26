using CSV, DataFrames, Plots, Plots.PlotMeasures, LaTeXStrings
cd(@__DIR__)
states = [1, 3, 5]
count = 0;
costf = collect(0.05:0.05:0.5)
T = collect(1:1:6);

#plotlyjs()
#pgfplotsx()
myplots = [];
for s=1:3
    sensitivity = Matrix(CSV.read("sensitivity.csv",DataFrame))[1+10*(s-1):10*s,:];

    p1 = plot3d(
        1, 
        xlim = (1,6),
        ylim=(0,0.5),
        zlim = (0.0,maximum(sensitivity)+5),
        legend=false,
        xtickfontsize = 17,
        ytickfontsize = 17,
        ztickfontsize = 17,
        size=(700,700),
        #guidefont=30,
        fmt = :png,
        #camera = (20,20),
        top_margin=-10mm,
        bottom_margin=-10mm,
        left_margin=-2mm,
        right_margin=-10mm,
        widen=true,
        gridlinewidth = 3,
        #xrotation=0,
        )
    clcount=1.0;
    clcount2=0.0;
    for j=1:10
        x = []; y = []; z = [];
        for i=1:6
            push!(x, T[i]); push!(y, costf[j]); push!(z, sensitivity[j,i])
        end
        #xlabel=L"t",
        #ylabel= L"\nu",
        #zlabel=L"\bar{f}",

        plot!(
            x,y,z,
            label = false,
            lc=RGB(clcount2,0,clcount),
            lw=7, 
            )

        #plot!([2.5], [0.01], text=L"t", ms=30)
        plot!([2.5], [0.03], text=text(L"t", 40))
        plot!([5.25], [0.055], text=text(L"\nu", 40))
        if s==1
            plot!([0.75], [0.225], text=text(L"\bar{f}", 40))
        else
            plot!([0.65], [0.225], text=text(L"\bar{f}", 40))
        end
        
        clcount -= 0.1
        clcount2 += 0.1
    end
    
    push!(myplots,p1)
    savefig(p1, "p"*string(states[s])*".png")
end
myplots[1]
