using Plots
using LaTeXStrings

function save(fig, savename)
    savefig(fig, savename*".png")
    savefig(fig, savename*".svg")
    savefig(fig, savename*".pdf")
end

function scatterplot(x1, x2,y,col,alpha; colorbar=true, xlabel="UMAP 1", ylabel="UMAP 2",clim=(-1,1))
    xmarks = [0, "", 0.2, "", 0.4, "", 0.6, "", 0.8, "", 1]
    xvalues=range(0,1, length=length(xmarks))
        
    if y==nothing
        pl = scatter(x1, x2, label="",c=col, alpha=alpha, clim=clim, markersize=8, markerstrokewidth=1, aspect_ratio=1, size=(450,450), colorbar=colorbar)
    else
        pl = scatter(x1, x2, label="",c=col, alpha=alpha, clim=clim, markersize=8, markerstrokewidth=1, aspect_ratio=1, size=(450,450), colorbar=colorbar,  marker_z=y)
    end
    plot!(xlim=(0,1),ylim=(0,1))
    plot!(guidefontsize=20, tickfontsize=15,right_margin = 12Plots.mm)
    plot!(xticks=(xvalues, xmarks), yticks=(xvalues, xmarks), xlabel=xlabel, ylabel =ylabel)
    return pl
end

function heatmapplot(matrix,col; colorbar=true, xlabel=L"x_1", ylabel=L"x_2",clim=(0,1))
    xmarks = [0, "", 0.2, "", 0.4, "", 0.6, "", 0.8, "", 1]
    xvalues=range(1,size(matrix,1), length=length(xmarks))

    pl = heatmap(matrix, c=col,clim=clim,colorbar=colorbar, aspect_ratio=:equal,size=(442,400), right_margin = 5Plots.mm, xlabel=xlabel, ylabel=ylabel)
    plot!(xticks=(xvalues, xmarks), yticks=(xvalues, xmarks))
    plot!(guidefontsize=20, tickfontsize=15,left_margin = 5Plots.mm)
    return pl
end

function upperlowerbound(pxy, y, α)
    l = []
    u = []
    C = cumsum(pxy, dims=2)
    C = (C.>α/2) .* (C.< 1-α/2)  #symmetric: same uncovered part left and right
    for i in 1:size(pxy,1)
        v = y[C[i,:].==1]
        if length(v)==0  #can still throw error, if C switches from 0 to 1 in one step. Very narrow prediction, not interval. In this case, should increase lenghtscale. Or increase D dimensionality.
            push!(l, y[argmax(pxy[i,:])])
            push!(u, y[argmax(pxy[i,:])])
        else
            push!(l, v[1])  
            push!(u, v[end])
        end
    end
    return l,u
end

function visualize_pxy(pxy, y, confidence=0.95)
    #denoise and normalize
    ϵ = abs(min(minimum(pxy), 0))
    pxy[abs.(pxy).<=ϵ].=0
    pxy = pxy ./ sum(pxy, dims=2)
    
    #compute properties
    EVE = pxy * y
    MLE = y[[argmax(pxy[i,:]) for i in 1:size(pxy,1)]]
    l,u = upperlowerbound(pxy, y, 1-confidence)
    
    xvalues = range(0,1,length=11);
    xmarks = [0, "", 0.2, "", 0.4, "", 0.6, "", 0.8, "", 1]

    pl = plot(xticks=(xvalues, xmarks), yticks=(xvalues, xmarks), xlabel=L"x", ylabel=L"y")
    plot!(x, Float64.(l), fillrange = Float64.(u), fillalpha = 0.1, c=:blue, linewidth=0, legend=false)
    plot!(x, ftrue, c=:black, label=L"$f$", linewidth=3)
    plot!(xs,fs, seriestype=:scatter, markersize=5, c=:green)
    plot!(x, MLE, c=:blue, linewidth=3,linestyle=:dashdotdot)
    plot!(x, EVE, c=:blue, linewidth=3,linestyle=:solid)
    plot!(guidefontsize=20, tickfontsize=15,left_margin = 5Plots.mm)
    return pl
end