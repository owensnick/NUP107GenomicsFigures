
maxnorm(X) = X./maximum(X, dims=2)

function kmeansseed(X, k, seed=1618)
    Random.seed!(seed)
    kmeans(X, k)
end

function kmeans_ud_silhouette(meta, tpm, filtind; ks=2:10, seed=1618)
    ui_labels = filter(n -> occursin("UC_", string(n)), names(tpm))
    mo_labels = filter(n -> occursin("MO_", string(n)), names(tpm))
    

    U = Matrix{Float64}(tpm[filtind, ui_labels])
    M = Matrix{Float64}(tpm[filtind, mo_labels])

    u_ge_m = vec(mean(U .- M, dims=2)) .≥ 0
    ind_u_ge_m = falses(length(filtind))
    ind_u_ge_m[filtind] .= u_ge_m

    
    UM_D = maxnorm([U M][  u_ge_m, :])
    UM_U = maxnorm([U M][.!u_ge_m, :])

    dist_D = pairwise(Euclidean(), UM_D')
    dist_U = pairwise(Euclidean(), UM_U')

    KMK_U =[kmeansseed(UM_U', k, seed) for k in ks]
    KMK_D =[kmeansseed(UM_D', k, seed) for k in ks]
    

    costs_u = [KMK.totalcost for KMK in KMK_U]
    costs_d = [KMK.totalcost for KMK in KMK_U]

    SD = @showprogress [silhouettes(KMK, dist_D) for KMK in KMK_D]
    SU = @showprogress [silhouettes(KMK, dist_U) for KMK in KMK_U]

    # SD = @showprogress [silhouettes(kmeans(UM_D', k), dist_D) for k in ks]
    # SU = @showprogress [silhouettes(kmeans(UM_U', k), dist_U) for k in ks]

    pu = boxplot(ks, SU, leg=false, ylabel="Silhouette", xlabel="k", legend=false, title="Up")
    iu = argmax(median.(SU))
    annotate!((ks[iu], 0.8, text("*")))

    pd = boxplot(ks, SD, leg=false, ylabel="Silhouette", xlabel="k", legend=false, title="Down")
    id = argmax(median.(SD))
    annotate!((ks[id], 0.8, text("*")))
    
    println("Optimal maximum K between U and D: k = ", ks[max(iu, id)])
    plot(pd, pu, xticks=ks, link=:y, size=(600, 300))
    
end


function kmeans_ud(meta, tpm, filtind, kd, ku; reorderU=Int[], reorderD=Int[], seed=1618)
    ui_labels = filter(n -> occursin("UC_", string(n)), names(tpm))
    mo_labels = filter(n -> occursin("MO_", string(n)), names(tpm))
    

    U = Matrix{Float64}(tpm[filtind, ui_labels])
    M = Matrix{Float64}(tpm[filtind, mo_labels])

    u_ge_m = vec(mean(U .- M, dims=2)) .≥ 0
    ind_u_ge_m = falses(length(filtind))
    ind_u_ge_m[filtind] .= u_ge_m


    UM_D = maxnorm([U M][  u_ge_m, :])
    UM_U = maxnorm([U M][.!u_ge_m, :])

    KMK_D = kmeansorder(UM_D', kd, seed)
    KMK_U = kmeansorder(UM_U', ku, seed)

    if !isempty(reorderU)
        KMK_U = reordercluster(KMK_U, reorderU);
    end
    if !isempty(reorderD)
        KMK_D = reordercluster(KMK_U, reorderD);
    end



    #### Make cluster table
    genename = last.(split.(tpm.Gene, "|"))
    gtable = DataFrame(Gene=tpm.Gene, GeneName=genename, Named=occursin.("|", tpm.Gene) .& .!occursin.(r"^LOC|^Xetrov", genename),
                       Ind=filtind, D_ind = ind_u_ge_m, U_ind = .!ind_u_ge_m, ClusterD=zeros(Int, length(filtind)), ClusterU=zeros(Int, length(filtind)))
    # @show size(gtable)
    # @show size(KMK_D.A)
    gtable.ClusterD[filtind .&   ind_u_ge_m] = KMK_D.A
    gtable.ClusterU[filtind .& .!ind_u_ge_m] = KMK_U.A
    gtable.Cluster = clabel.(gtable.ClusterD, gtable.ClusterU)

    addtable = tpm[!, Not(r"^Gene|^UIC_|^hiK_|^CR_")]
    if !isempty(addtable)
        gtable = [gtable addtable]
    end

    
    (KMK_D=KMK_D, KMK_U=KMK_U, gtable=gtable, meta=meta)

end

function clabel(cd, cu)
   
    if (cd == 0) && (cu == 0)
        return "N"
    elseif (cd == 0) && (cu != 0)
        return string("U", cu)
    elseif (cd != 0) && (cu == 0)
        return string("D", cd)
    else
       error("CL Assignment: $cd, $cu") 
    end
        
end

function clusterorder(cl)
    if cl == "N"
        return (0, '0')
    elseif cl[1] == 'U'
        return (1, cl[2])
    elseif cl[1] == 'D'
        return (2, cl[2])
    else
       error("CL unknown: $cl") 
    end
end


function join_cluster_tables(KMK_UD_PA, KMK_UD_RZ)
    @assert KMK_UD_PA.gtable.Gene == KMK_UD_RZ.gtable.Gene
    gt = KMK_UD_PA.gtable[!, [:Gene, :GeneName, :Named]]
    pa = KMK_UD_PA.gtable[!, [:Ind, :D_ind, :U_ind, :ClusterD, :ClusterU]]
    rz = KMK_UD_RZ.gtable[!, [:Ind, :D_ind, :U_ind, :ClusterD, :ClusterU]]
    
    
    pa.Cluster = clabel.(pa.ClusterD, pa.ClusterU)
    rz.Cluster = clabel.(rz.ClusterD, rz.ClusterU)
    rename!(pa, string.("PA_", names(pa)))
    rename!(rz, string.("RZ_", names(rz)))
    
    
    [gt pa rz]
end


##### Cluster plotting functions

function plotclusterud_stack(KMK_UD; lta = @layout [a ; b{0.35h}])


    phu = viscluster_double(KMK_UD.meta.Time[KMK_UD.meta.Treat .== "UC"], KMK_UD.KMK_U, c=:viridis, colorbar=false, ylabel="Up", titlefont=font(9, "helvetica"), xticks=false, bottom_margin=-2mm, title="UC               MO")
    phd = viscluster_double(KMK_UD.meta.Time[KMK_UD.meta.Treat .== "UC"], KMK_UD.KMK_D, c=:viridis, colorbar=false, ylabel="Down", titlefont=font(9, "helvetica"), top_margin=-1mm, xlabel="Time (hpf)", bottom_margin=3mm)
    plot!(phu, xticks=false)
    pcd = plotclustermeansd(KMK_UD.KMK_D, KMK_UD.meta, titlelabel="Down", xlabel="Time (hpf)", bottom_margin=3mm, minplot=3)
    pcu = plotclustermeansd(KMK_UD.KMK_U, KMK_UD.meta, titlelabel="Up",minplot=3)
    plot!(pcd, xticks=2.5:2.5:12.5)
    plot!(pcu, xticks=2.5:2.5:12.5)

    lt = @layout  [[a{0.3w} ; b{0.3w, 0.3h}] [c ; d]]
    # lt = @layout  [a{0.3w} b ; c{0.3h} d]
    
    ltb = @layout [a{0.2w} b{0.05w} [ c; d]]
    pa = plot(phu, phd, layout=lta)

    h2 = scatter([0,0], [0,1], zcolor=[0,1], clims=(0, 1),  xlims=(1,1.1), label="", c=:viridis, framestyle=:none)
    p = plot(pa, h2, pcu, pcd, layout=ltb, size=(900, 300), xguidefont=font(8, "helvetica"))
end


function viscluster_double(xp, KMK ; dy=2, kwargs...)
    X = KMK.X'
    delta = xp[2] - xp[1]
    dxp = [xp ; xp[end] .+ delta .+ xp]
    
    heatmap(dxp, 1:dy:size(X, 1), average_heatmap(X[KMK.SI, :], 1, dy), yflip=true; kwargs...)
    
    vline!([last(xp) .+ delta .+ 1.0], c=:white, lab="")
    clusterlines!(KMK)

    tp = [5.0, 10.0, xp[end] + delta + 5.0, xp[end] + delta + 10.0]
    
    xticks!(tp, ["5", "10", "5", "10"])


end

function clusterlines!(KMK)
    cn = KMK.KM.counts[KMK.KSI]
    hline!(cumsum(cn), c=:white, lab="")
    yticks!(cumsum(cn) .- cn/2, string.(1:KMK.k))
end


function plotclustermeansd(KMK, meta; X = KMK.X, titlelabel = "", labels = string.(titlelabel, meta.Label), lt = :row2, clusterby=["U", "M"], timescale=true, stageticks=true, minplot=4, kwargs...)

    # pointwise mean and sd M and SD
    phs = Plots.Plot[]
    cc = [:steelblue, :orange]
    for k = 1:KMK.k
        
        M = vec(mean(X[:, KMK.A .== k], dims=2))
        S = vec( std(X[:, KMK.A .== k], dims=2))
    
        p = plot(grid=false, title=string(titlelabel, " C", k, " : ", sum(KMK.A .== k)), titlefont=font(9, "helvetica"), framestyle=:box; kwargs...)
        for (i, cl) in enumerate(clusterby)
            if cl == "U"
                ind = meta.Treat .== "UC"
            elseif cl == "M"
                ind = meta.Treat .== "MO"
            end
            
            
            t = meta.Time[ind]
            
            plot!(t, M[ind], lab=ifelse(k == 1, cl, ""), ribbon=S[ind],  fillcolor=cc[i], c=cc[i], xlims=extrema(t))
        end
        
        push!(phs, p)
    end

    if length(phs) < minplot
        push!(phs, plot(axis=false, framestyle=nothing, grid=false, xticks=false))
    end

    [plot!(p, left_margin=-1mm, right_margin=-1mm, yformatter=y->"", yticks=false) for p in phs[2:end]]
    plot!(phs[1], right_margin=-1mm, yticks=0:1)

    plot(phs..., layout=(1, length(phs)), ylims=(0, 1.05))
end

function annotate_parz_clusterplot!(p = plot!(); kwargs...)
    vline!([1.5], c=:white, ls=:dash, lab="")
    hline!([1.5], c=:white, ls=:dash, lab="")
    vline!([4.5], c=:white, ls=:dash, lab="")
    hline!([4.5], c=:white, ls=:dash, lab="")
    xlabel!("PA Cluster")
    ylabel!("RZ Cluster")
end


function plotclustermeansd_vert(KMK, meta; X = KMK.X, titlelabel = "", labels = string.(titlelabel, meta.Label), lt = :row2, clusterby=["U", "M"], timescale=true, stageticks=true, minplot=3, kwargs...)

    # pointwise mean and sd M and SD
    phs = Plots.Plot[]
    cc = [:steelblue, :orange]
    for k = 1:KMK.k
        
        M = vec(mean(X[:, KMK.A .== k], dims=2))
        S = vec( std(X[:, KMK.A .== k], dims=2))
    
        p = plot(grid=false, ylabel=string(titlelabel, " C", k, " : ", sum(KMK.A .== k)), titlefont=font(9, "helvetica"), yticks=0:1, ylims=(0, 1.05), framestyle=:box; kwargs...)
        for (i, cl) in enumerate(clusterby)
            if cl == "U"
                ind = meta.Treat .== "UC"
            elseif cl == "M"
                ind = meta.Treat .== "MO"
            end
                
            t = meta.Time[ind]
            xts = 1:2:13
            if timescale && !stageticks
                t = t./2 .+ 4
                xts = 4:2:10
            end
        

            plot!(t, M[ind], lab=ifelse(k == 1, cl, ""), ribbon=S[ind],  fillcolor=cc[i], c=cc[i], xlims=extrema(t))
        end
        
        push!(phs, p)
    end

    if length(phs) < minplot
        push!(phs, plot(axis=false, framestyle=nothing, grid=false, xticks=false))

    end

    plot(phs..., layout=(length(phs), 1)), phs

end