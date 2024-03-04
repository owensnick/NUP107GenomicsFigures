function savedisplay(f; dir="", projdir=getprojectdir(), kwargs...)
    path = joinpath(projdir, "figs", dir)
    if !isdir(path)
        mkpath(path)
    end
    p = plot!(fontfamily="helvetica")
    display(p)
    file = joinpath(path, f)
    savefig(file; kwargs...)
end

tt(x) = string(x)
tt(x::Float64; digits=1) = string(round(x, digits=digits))


function heatannot(xlabs, ylabs, M; trf=identity, fs=font("helvetica", 9), annot=true, kwargs...)
    p = heatmap(1:size(M, 2), 1:size(M, 1), trf.(M), xticks=(1:size(M, 2), xlabs), yticks=(1:size(M, 1), ylabs), yflip=false; kwargs...) 
    annot && annotate!([(i, j, text(tt(M[j, i]), fs)) for i = 1:size(M, 2), j = 1:size(M, 1)][:])
    p
end

function plotnupgene(gene, label, nups; tmax = 24, ps=1e+6, p = plot(fontfamily="helvetica"), kwargs...)
    nsel = @subset(nups, occursin.(gene, :Gene), :Label .== label, :Time .<= tmax)
    
    data = @subset(nsel, :DataType .== "data")
    med  = @subset(nsel, :DataType .== "gpmed")
    ciu  = @subset(nsel, :DataType .== "gpciu")
    cil  = @subset(nsel, :DataType .== "gpcil")
    
    @assert issorted(med.Time)
    @assert med.Time == ciu.Time
    @assert med.Time == cil.Time
    
    plot!(p,  kwargs...)
    scatter!(data.Time, data.TPE/ps, marker=:circle, markersize=1, c=:black, grid=false, leg=false)

    plot!(med.Time, med.TPE/ps, ribbon=[med.TPE .- cil.TPE ciu.TPE .- med.TPE]./ps, c=:black, fillalpha=0.2)
    
    
    yl = ylims()
    ylims!(0, yl[2])
    
    annotate!(10, yl[2], text(gene, font("helvetica", 10)))
    plot!(yticks=ytm(yl[2]))

    
end

function ytm(ym)
    (ym > 8) && return [0.0, 8.0]
    (ym > 5) && return [0.0, 5.0]
    (ym > 2) && return [0.0, 2.0]
    (ym > 1) && return [0.0, 1.0]
    return [0.0, 0.5]
   
end

function nupgeneplot(nups, parz, tmax=24, xt=10:10:20)
    nsel = @subset(nups, :Label .== parz, :Time .<= 24)
    med = @subset(nsel, :DataType .== "gpmed")
    st = sort(unique(med.Time))
    
    M = unstack(med, :Gene, :Time, :TPE)
    dropmissing!(M)
     
    
    X = Matrix(M[!, string.(st)])
    X ./= maximum(X, dims=2)
    KMK = kmeansorder(X', 6)

    genes = last.(split.(M.Gene[KMK.SI], "|"))

    ph = heatmap(st, 1:size(X, 1), X[KMK.SI, :], yticks=(1:size(X, 1), genes), fontfamily="helvetica")
    phs = [plotnupgene(g, parz, nsel, tmax=tmax) for g in genes]
    p = plot(phs..., size=(800, 500), leg=false, leftmargin=-2mm, rightmargin=-1mm, bottommargin=-1mm, titlefont=font(10, "helvetica"), xticks=xt)
    plot(ph, p, layout=grid(1, 2, widths=[0.3, 0.7]), size=(800, 400)) 
end


function mir_enrich_heatmap(me_all; pt = 0.01, fs=font("helvetica", 9), kwargs...)
    
    us = unstack(me_all, [:RCSeed, :LevenshteinMir427, :IDs], :Cluster, :pvr)
    us = us[!, [names(us)[1:3] ; sort(names(us)[4:end])]]
    
    if size(us, 1) < 4
       return plot() 
    end
    
    P = -log10.(coalesce.(Matrix(us[!, 4:end]), 1))
    v = vec(maximum(P, dims=2)) .> -log10.(pt)
    if sum(v) > 4
        Z = P[v, :]
        # KMK = kmeansorder(Z', 4)
        # labels = @with us string.(replace.(joincomma.(sort.(split.(:IDs, ","))), r"xtr-|miR-" => ""), " ", :LevenshteinMir427, " ", :RCSeed)
        labels = @with us string.(replace.(:IDs, r"xtr-|miR-" => ""), " ", :LevenshteinMir427, " ", :RCSeed)
        sdf = DataFrame(U=us.LevenshteinMir427[v], Z=Z[:, 4])#, KSI=KMK.SI)
        si = sortperm(sdf)
        cs = combine(groupby(sdf, :U, sort=true), nrow => :count)
        
        heatannot(names(us)[4:end], labels[v][si, :], abs.(Z[si, :]), annot=true; yflip=true, fontfamily="helvetica", fs=fs, kwargs...)
        hline!(cumsum(cs.count) .+ 0.5, c=:white, lab="", ls=:dash)
    else
       plot() 
    end
end

function utr_motif_heatmap(moten; kwargs...)
    us = unstack(moten, [:ID], :Cluster, :pvr)
    us = us[!, ["ID" ;  sort(names(us)[2:end]) ]]
    P = -log10.(coalesce.(Matrix(us[!, 2:end]), 0)) .+ 0

    heatannot( names(us)[2:end], us.ID, P; kwargs...)
end


plotgpset(args...; kwargs...) = (plot(); plotgpset!(args...; kwargs...))
function plotgpset!(gene, genes, gpr, altids=[]; plotall=false, layout=[], kwargs...)
        
    ind = findall(occursin.(gene, genes))
    
    isempty(ind) && error("$gene not found")
    if length(ind) > 1
        println("Multiple matches for $gene:\n$(genes[ind])")
        phs = [plotgpset(gpr[i]; suptitle=genes[i], kwargs...) for i in ind]
        if isempty(layout)
            plot(phs..., size=(1200, 500); kwargs...)
        else
            plot(phs..., layout=layout, size=(1200, 500); kwargs...)
        end
        
    else
        plotgpset!(gpr[first(ind)] ; suptitle=genes[first(ind)], kwargs...)
    end
end



function plotgpset!(gpr; samples=["U", "M"], showparams=false, suptitle="", kwargs...)
    p = plot(; kwargs...)
    cc = Dict("U" => :steelblue, "M" => :orange, "F" => :black)
    for s in samples
        gp = getindex(gpr, Symbol(lowercase(s)))
        plotgp!(gp, c=cc[s], lab=s, showparams=showparams; kwargs...)
    end
    

    if showparams
        lml_strs = String[]
        param_strs = String[]
        for s in samples
            gp = getindex(gpr, Symbol(lowercase(s)))

            push!(lml_strs, string("lml(", s, ") = ", tt(gp.gp.mll)))
            push!(param_strs, gpparamstring(gpr, Symbol(lowercase(s))))
        end
        
        s = string(" ", join(lml_strs, " | "), "\n",
                   "lr, bic = ", tt.(gp_bic(gpr.f, [gpr.u, gpr.m])), "\n",
                   join(param_strs, "\n"), "\n")
        plot!(titlefont=font(8), top_margin=10mm)
    else
        lr, bic = gp_bic(gpr.f, [gpr.u, gpr.m])
        md = maxcohensd(gpr)
        s = string("  LR: ", tt(lr), " CD: ", tt(md))
        plot!(titlefont=font(10, "helvetica"))
    end
    yl = ylims()
    ylims!(0, yl[2])
    plot!(title=string(suptitle, s))
    
end

plotgp(args...; kwargs...) = (plot(); plotgp!(args...; kwargs...))
tt(x) = string(x)
tt(x::Float64; digits=1) = string(round(x, digits=digits))
function gpparamstring(gpr, field)
    string("p(", field, ") = (",
     tt(gpr[field].gp.kernel.σ2), ", ",
     tt(gpr[field].gp.kernel.ℓ), ", ",
     tt(exp(2*gpr[field].gp.logNoise.value)), ")")
end

function plotgp!(gpr; c=:auto, marker=:auto, plotribbon=true, showparams=false, lab="", kwargs...)


    ### to plot gpr.f assume that gpr.t is arranged [t ; t]
    n = length(gpr.t) ÷ 2
    if gpr.t[n+1] > gpr.t[n]
        plot!(gpr.t, gpr.y, marker=(:circle), c=c, lab=""; kwargs..., alpha=0.75)
    else
        scatter!(gpr.t[1:n], gpr.y[1:n], marker=(:steelblue, 2), c=c, lab="U data"; kwargs...)
        scatter!(gpr.t[(n+1):end], gpr.y[(n+1):end], marker=(:orange), 2, c=c, lab="M data"; kwargs...)
    end

     if plotribbon
         p = plot!(gpr.st, gpr.sf, ribbon=[gpr.sf - gpr.cil gpr.ciu - gpr.sf], c=c, lab=lab; kwargs...)
     else
         p = plot!(gpr.st, gpr.sf,  c=c, lab=lab, kwargs...)
     end

     if showparams
        s = string("LML = ", tt(gpr.gp.mll), "\np = (",
         tt(gpr.gp.kernel.σ2), ", ",
         tt(gpr.gp.kernel.ℓ), ", ",
         tt(exp(2*gpr.gp.logNoise.value)), ")")
        plot!(title=s, titlefont=font(8, :white))
     end
     p
end