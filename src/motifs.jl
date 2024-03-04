function loadmotifselection(; projdir=getprojectdir())
    path = joinpath(projdir, "data", "motifs")
    files = glob("*.meme", path)
    motsel = loadmeme.(files, 1/100, Ref(MotifScanner.background_xenopus_zero_markov()))
    motsel
end



function loadpromoters(file="XENTR_9.1_Xenbase_spike.nmt.prom500.tsv.gz", idd=genenamedict(), projdir=getprojectdir())
    filepath = joinpath(projdir, "data", file)
    proms = CSV.read(filepath, DataFrame)


    gidns = get_gene_id_name.(proms.GeneID, Ref(idd))
    ids = first.(gidns)
    genenames = last.(gidns)
    genes = string.(ids, "|", genenames)
    proms.Gene = genes
    proms.Isoform = string.(ids, "|", proms.TranscriptID, "|", genenames)
    proms.Index = 1:size(proms, 1)


    proms = proms[!, [:Gene, :Isoform, :GeneName, :GeneID, :TranscriptID, :strand, :Index, :chrom, :start, :stop, :TSS, :PromoterStart, :PromoterStop, :Length, :seq]]
    proms
    proms = @subset(proms, :Gene .∈ Ref(Set(gex_pa_rz.gpgenes)), :Length .== 500)
    proms.Index = 1:size(proms, 1);
    proms.seq = replace.(proms.seq, "N" => "-")
    proms.seq = LongDNA{4}.(uppercase.(proms.seq))
    gene_proms = combine(groupby(proms, [:Gene, :strand]), :Index => Ref => :Index);

    proms, gene_proms
end


function motenrightment_clusters()
    clsel = [["U1"], ["U2"], ["U3"], 
    ["D1"], ["D2"], ["D3"],
    ["U3", "D3"],
    ["U1", "U2", "U3", "D1", "D2", "D3"]]
    clsel
end

function scanmotifmax(mots, seqs)
    M = Matrix{Float64}(undef, length(seqs), length(mots))
    
    n = length(mots)
    p = Progress(n, .5, "Scanning Motifs")
    Threads.@threads for i = 1:n
        for j = 1:length(seqs)
            M[j, i] = MotifScanner.scanmax(seqs[j], mots[i]) |> first
        end
        next!(p)
    end
    M
end

function mot_score_thresh_enrich(mots, proms, gene_proms, parz_clusters; uc=motenrightment_clusters(), motlen=100, minscore=5)
    isoscores =  scanmotifmax(mots, proms.seq)
    scores = collapse_isoscores(gene_proms.Index, isoscores, collapse_fun=maximum)
    labels = getindex.(mots, :id)
    scoredf = coalesce.(leftjoin(parz_clusters, [DataFrame(Gene=gene_proms.Gene) DataFrame(scores, labels)], on=:Gene), 0)
    score_thresh = Vector{Float64}(undef, length(mots))    
    res = Vector{Base.return_types(ht, (BitVector, BitVector))[1]}()
    @showprogress for (k, (m, l)) in enumerate(zip(mots, labels))
        mt = range(minscore, maximum(m.pbg, dims=1) |> sum, length=motlen)
        htr = [ht(scoredf.PA_Cluster .∈ Ref(Set(u)), scoredf[!, l] .> sc) for u in uc, sc in mt]
        pvr = getindex.(htr, :pvr)
        i, j = argmin(pvr).I
        append!(res, htr[:, j])
        score_thresh[k] = mt[j]
            
    end
    scoredf = DataFrame(MotifID=labels, MotifName=getindex.(mots, :name), ScoreThresh=score_thresh)
    ht_df = [DataFrame(MotifID=repeat(labels, inner=length(uc)), MotifName=repeat(getindex.(mots, :name), inner=length(uc)), ClusterGroup=repeat(uc, length(labels))) DataFrame(res)]
    
    innerjoin(scoredf, ht_df, on=[:MotifID, :MotifName])
    
end

