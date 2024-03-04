
function loadnupgenes(;file="xep_parz_nupgenes_tpe.tsv.gz", projdir=getprojectdir())
    filepath = joinpath(projdir, "data", file)
    nups = CSV.read(filepath, DataFrame)
    nups
end


function genenamedict(; file="XENTR_9.1_Xenbase_spike.nmt.ids.tsv.gz", projdir=getprojectdir())

    filepath = joinpath(projdir, "data", file)
    ids = CSV.read(filepath, DataFrame)

    idd = Dict{String, String}()

    for (id, gn) in zip(ids.GeneID, ids.GeneName)
        (id == gn) && continue
        if haskey(idd, id) 
            if idd[id] != gn
                error("Multiple matches for $id:\tgn=$gn\t$(idd[id])")
            end
        else
            idd[id] = gn
        end

    end
    idd

end

# function load_trans_ids(; file="XENTR_9.1_Xenbase_spike.nmt.ids.tsv", projdir=getprojectdir())

#     filepath = joinpath(projdir, "data", file)
#     ids = CSV.read(filepath, DataFrame)

#     iddict = Dict{String, String}()

#     for (trid, gid, gn) in zip(ids.TranscriptID, ids.GeneID, ids.GeneName)
#         gene = string(gid, "|", gn)
#         if haskey(iddict, trid)
#             if gene != iddict[trid]
#                 error("Multiple matches for $trid: $gene $(iddict[trid])")
#             end
#         end
#         iddict[trid] = gene
#     end
#     iddict
# end


function load_meta_rsem(;rsemdir="aligns-rsemstar-xt9", projdir=getprojectdir())
    pth = joinpath(projdir, "data", rsemdir )
    files = glob("*.genes.results", pth)
    fields = split.(basename.(files), "_")
    treat = ifelse.(startswith.(getindex.(fields, 2), "UIC"), "UC", "MO")
    time = parse.(Float64, replace.(getindex.(fields, 2), Ref(r"UIC|MO" => ""))) .+ ifelse.(getindex.(fields, 3) .== "5", .5, .0)
    labels = string.(treat, "_", replace.(string.(round.(time, digits=1)), Ref("." => "_")))
    meta = DataFrame(Study="Nup107", Treat=treat, Time=time, Label=labels, File = files)
    meta = sort(meta, [order(:Treat, by=x -> Dict("UC" => 1, "MO" => 2)[x]), :Time])
    meta[!, :Index] = 1:size(meta, 1)
    meta
end



function load_rsem_isoforms(meta, idd=genenamedict(); field = :TPM)
    isofiles = replace.(meta.File, ".genes." => ".isoforms.")
    @assert all(isfile, isofiles)
    tables = @showprogress [CSV.read(f, DataFrame) for f in isofiles]
    geneids = first(tables).gene_id
    isoids =  first(tables).transcript_id

    @assert all(df -> df.gene_id == geneids, tables)
    @assert all(df -> df.transcript_id == isoids, tables)

    # get gene id from isofield
    isofields = split.(isoids, "|")
    geneids = first.(isofields)
    
    gidns = get_gene_id_name.(geneids, Ref(idd))
    ids = first.(gidns)
    genenames = last.(gidns)
    genes = string.(ids, "|", genenames)
    
    visoids = string.(ids, "|", getindex.(isofields, 2), "|", genenames)


    isotpm = [DataFrame(Gene=genes, Isoform=visoids) DataFrame(getproperty.(tables, field), meta.Label)]
    isotpm

end

function get_gene_id_name(gene, idd=genenamedict())
    if occursin("_", gene)
        fields = split(gene, "_")
        geneid = fields[1]
        genename = fields[2]
        
        if haskey(idd, geneid)
            if genename != idd[geneid]
                # println("Multiple ids: ", gene, "\t", idd[geneid])
            end
        end
        geneid, genename
    else
        gene, get(idd, gene, "")
    end
end

function load_rsem_tpm(meta, idd=genenamedict(); field=:TPM)

    tables = @showprogress [CSV.read(f, DataFrame) for f in meta.File]
    geneids = first(tables).gene_id
    @assert all(df -> df.gene_id == geneids, tables)

    gidns = get_gene_id_name.(geneids, Ref(idd))
    ids = first.(gidns)
    genenames = last.(gidns)
    
    ### check that all genenames are unique
    # genids = @subset(DataFrame(Gene=ids, GeneName=genenames), :GeneName .!= "")
    # multnames = @subset(combine(groupby(genids, :GeneName), nrow => :Count, :Gene => Ref => :GeneIDs), :Count .> 1)
    # multnames |> display
    # @with multnames println.(:GeneName, "\t", :GeneIDs);
    # non_blank_names = filter(!isempty, genenames)

    # @assert length(unique(non_blank_names)) == length(non_blank_names) # no duplicates

    genes = string.(ids, "|", genenames)
    # genes = string.(geneids, "|", get.(Ref(idd), geneids, ""))

    tpm = [DataFrame(Gene=genes) DataFrame(getproperty.(tables, field), meta.Label)]
    tpm
end




function calctotalreads(meta)
    totals = @showprogress [sum(CSV.read(f, DataFrame).expected_count) for f in meta.File]
    DataFrame(Label=meta.Label, TotalReads=totals)
end

function longest_run(x, τ = 0)
    rl = 0
    mrl  = 0
    @inbounds for v ∈ x
        rl  = ifelse(v > τ, rl + 1, 0)
        mrl = max(mrl, rl)
    end
    mrl
end


function tablestats_filter(meta, tpm; lrt=0.4, rl = 6, ns=2)

    T = Matrix(tpm[!, meta.Label])

    minv  = vec(minimum(T, dims=2))
    maxv  = vec(maximum(T, dims=2))
    meanv = vec(mean(T, dims=2))
    minv  = vec(minimum(T, dims=2))
    
    ui_ind = meta.Treat .== "UC"
    mo_ind = meta.Treat .== "MO"


    lr_UI = [longest_run(c[ui_ind], lrt) for c in eachrow(T)]
    lr_MO = [longest_run(c[mo_ind], lrt) for c in eachrow(T)]

    stats = DataFrame(Gene=tpm.Gene, minv=minv, meanv=meanv, maxv=maxv, lr_UI=lr_UI, lr_MO=lr_MO, Index=1:length(minv))

    filtind = mapreduce(l -> stats[!, l] .>= rl, +, [:lr_UI, :lr_MO]) .>= ns


    stats, filtind
end




function load_tpm_pa_rz(padir="aligns-rsemstar-xt9", rzdir="ribozero-aligns-rsemstar-xt9")
    meta_pa = load_meta_rsem(rsemdir=padir);
    tpm_pa = load_rsem_tpm(meta_pa);
    stats_pa, filtind_pa = tablestats_filter(meta_pa, tpm_pa);
    
    exc_pa = load_rsem_tpm(meta_pa, field=:expected_count)
    total_counts_pa = describe(exc_pa[!, meta_pa.Label], sum => :TotalCounts)
    total_counts_pa.variable = string.(total_counts_pa.variable)
    @assert meta_pa.Label == total_counts_pa.variable
    meta_pa.TotalCounts = total_counts_pa.TotalCounts

    
    normcounts_pa = norm_exc_counts(exc_pa)
    stats_nc_pa, filtind_nc_pa = tablestats_filter(meta_pa, normcounts_pa, lrt=10);
    
    meta_rz = load_meta_rsem(rsemdir=rzdir);
    tpm_rz = load_rsem_tpm(meta_rz);
    stats_rz, filtind_rz = tablestats_filter(meta_rz, tpm_rz);
    
    exc_rz = load_rsem_tpm(meta_rz, field=:expected_count)
    total_counts_rz = describe(exc_rz[!, 2:end], sum => :TotalCounts)
    total_counts_rz.variable = string.(total_counts_rz.variable)
    @assert meta_rz.Label == total_counts_rz.variable
    meta_rz.TotalCounts = total_counts_rz.TotalCounts
    # meta_rz = innerjoin(meta_rz, total_counts_rz, on=:Label => :variable)
    # sort!(meta_rz, :Index)
    normcounts_rz = norm_exc_counts(exc_rz)
    stats_nc_rz, filtind_nc_rz = tablestats_filter(meta_rz, normcounts_rz, lrt=10);
    
    ### ensure agreement and compare filtering
    @assert tpm_rz.Gene == tpm_pa.Gene
    @assert meta_rz.Time == meta_pa.Time
    @assert meta_rz.Treat == meta_pa.Treat

    
    combine(groupby(DataFrame(Data="TPM", FIPA=filtind_pa, FIRZ=filtind_rz), [:FIPA, :FIRZ]), nrow => :count) |> display
    combine(groupby(DataFrame(Data="NormCounts", FIPA=filtind_nc_pa, FIRZ=filtind_nc_rz), [:FIPA, :FIRZ]), nrow => :count) |> display

    filtind = filtind_pa .& filtind_rz
    f_ind = findall(filtind);
    gpgenes = tpm_pa.Gene[filtind];
    
    filtind_nc = filtind_nc_pa .& filtind_nc_rz
    f_ind_nc = findall(filtind_nc);
    gpgenes_nc = tpm_pa.Gene[filtind_nc];

    
    
    (meta_pa=meta_pa, tpm_pa=tpm_pa, meta_rz=meta_rz, tpm_rz=tpm_rz, filtind=filtind, f_ind=f_ind, gpgenes=gpgenes, normcounts_pa=normcounts_pa, normcounts_rz=normcounts_rz, filtind_nc=filtind_nc, f_ind_nc=f_ind_nc, gpgenes_nc=gpgenes_nc)
end



function pa_rz_cor_heatmap(gex; meta=gex.meta_pa, field=:tpm, dt=2, kwargs...)
    
    filtind = ifelse(field == :tpm, gex.filtind, gex.filtind_nc)
    
    pa_field = Symbol(string(field, "_pa"))
    rz_field = Symbol(string(field, "_rz"))
    
    tpm_pa = getindex(gex, pa_field)
    tpm_rz = getindex(gex, rz_field)
    
    TPA = Matrix(tpm_pa[filtind, gex.meta_pa.Label]);
    TRZ = Matrix(tpm_rz[filtind, gex.meta_rz.Label]);
    
    
    n = size(gex.meta_pa, 1)
    ns = 20
    xylab = string("PA:UIC", " "^ns, "PA:MO", " "^ns, "RZ:UIC", " "^ns, "RZ:MO")
    heatmap(corspearman([TPA TRZ]), ticks=(1:dt:(2*n), string.([meta.Time[1:dt:n] ; meta.Time[1:dt:n]])), xlabel=xylab, ylabel=xylab, size=(800, 600); kwargs...)
    vline!([size(TPA, 2)] .+ 0.5, c=:white, lab="", lw=2, ls=:dash)
    hline!([size(TPA, 2)] .+ 0.5, c=:white, lab="", lw=2, ls=:dash)
    
end

function norm_exc_counts(exc; genefields=1, scale=20e+6) #(mean(totalreads_pa), mean(totalreads_rz)) = (22.539569425952383, 17.593786890238096)
    normcounts = [exc[!, 1:genefields] DataFrame([scale*c/sum(c) for c in eachcol(exc[!, (genefields+1):end])], names(exc)[(genefields+1):end])]
    
    normcounts
end


function pcaplot(gex; parz="pa", dir=1, kwargs...)
    meta = getindex(gex, Symbol("meta_", parz))
    tpm = getindex(gex, Symbol("tpm_", parz))

    T = Matrix(tpm[gex.filtind, meta.Label])

    pcaplot(meta.Time, meta.Treat, T, pcA=1, pcB=2, title=string("PCA ", parz), dir=dir; kwargs...)


end

function pcaplot(labels, group, X; pcA = 1, pcB = 2, dir=1, kwargs...)
    LX = log10.(X .+ 1/10);
    pcaD = fit(PCA, LX)
    TD = MultivariateStats.transform(pcaD, LX)'
    explained = principalvars(pcaD)/tvar(pcaD)

    
    scatter(dir*TD[:, 1], TD[:, 2], zcolor=labels, marker=:auto, group=group, leg=:outertopright,
    xlabel="PC $pcA: $(round(100*explained[pcA], digits=2)) %",
    ylabel="PC $pcB: $(round(100*explained[pcB], digits=2)) %")
    annotate!([(dir*TD[i, 1] + 10, TD[i, 2], text(labels[i], font(8))) for i = 1:21])
    plot!(grid=false; kwargs...)
end


function writetpm(gex_pa_rz; projdir=getprojectdir())
    path = joinpath(projdir, "results")
    mkpath(path)
    CSV.write(joinpath(path, "nup107_mo_tpm_polyA.tsv.gz"), gex_pa_rz.tpm_pa, delim='\t', compress=true)
    CSV.write(joinpath(path, "nup107_mo_tpm_ribozero.tsv.gz"), gex_pa_rz.tpm_pa, delim='\t', compress=true);
end