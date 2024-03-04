
levenshtein(seed::String, seq::String) = levenshtein(DNAKmer{length(seed)}(LongDNA{4}(seed)), LongDNA{4}(seq))
function levenshtein(seed::T, seq) where T
    d = length(seed)
    
    for (i, k) in EveryKmer{T}(seq)
        d = min(d, mismatches(k, seed))
    end
    d
end

joincomma(s) = join(s, ", ")

function loadmirs(file="mature.fa.gz"; projdir=getprojectdir(), organism=r"Homo|Xenopus", mir427_rc_seed = "AAAGCACTTTC")
    mirrecords = collect(FASTA.Reader(GzipDecompressorStream(open(joinpath(projdir, "data", file)))))

    mirs = DataFrame(ID=identifier.(mirrecords), Description=description.(mirrecords), Seq=replace.(FASTX.sequence.(mirrecords), "U" => "T"))
    mirs.Organism = [join(split(d)[3:end-1], " ") for d in mirs.Description]

    mirs = @subset(mirs, occursin.(organism, :Organism))
    mirs.Num = getindex.(split.(mirs.ID, "-"), 3)
    mirs.Length = length.(mirs.Seq)
    sort!(combine(groupby(mirs, :Organism), nrow => :count), :count)
    
    mirs.Seed = [s[2:7] for s in mirs.Seq]
    mirs.RCSeed = string.(reverse_complement.(LongDNA{4}.(mirs.Seed)))
    
    mirs.Seed7 = [s[2:8] for s in mirs.Seq]
    mirs.RCSeed7 = string.(reverse_complement.(LongDNA{4}.(mirs.Seed7)))

    mirs.LevenshteinMir427 = levenshtein.(mirs.RCSeed, mir427_rc_seed)


    mirgroups = combine(groupby(mirs, [:RCSeed, :LevenshteinMir427]), nrow => :count, :ID => joincomma => :IDs, :Organism => joincomma => :Organisms, :Num => joincomma => :Nums)
    mirgroups.MirIndex = 1:size(mirgroups, 1)
    mirs, mirgroups
end



function loadutrs(file="XENTR_9.1_Xenbase_spike.nmt.utrs.tsv.gz"; gex_tested = [], idd=genenamedict(), projdir=getprojectdir())
    filepath = joinpath(projdir, "data", file)
    utrs = CSV.read(filepath, DataFrame)


    gidns = get_gene_id_name.(utrs.GeneID, Ref(idd))
    ids = first.(gidns)
    genenames = last.(gidns)
    genes = string.(ids, "|", genenames)
    utrs.Gene = genes
    utrs.Isoform = string.(ids, "|", utrs.TranscriptID, "|", genenames)
    utrs.Index = 1:size(utrs, 1)
    utrs = utrs[!, [:Gene, :Isoform, :GeneName, :GeneID, :TranscriptID, :strand, :Index, :Length, :seq]]

    if !isempty(gex_tested)
        genes_tested = intersect((g.gpgenes for g in gex_tested)...)
        utrs = @subset(utrs, :Gene .∈ Ref(Set(genes_tested)))
    end
    utrs.Index = 1:size(utrs, 1)
    gene_utrs = combine(groupby(utrs, [:Gene, :strand]), :Index => Ref => :Index);
    utrs, gene_utrs
end



function load_motifs_of_interest()
    sois_dict = Dict("EDEN" => "UAUAUAUGUGUGUCUAUCGUCACUUGUAUGUCAAAUAUU", "ARE" => "AUUUA", "YTHDF2" => "RRACH")
    ks = sort(collect(keys(sois_dict)))
    sois = DataFrame(ID=ks, Seq=replace.(getindex.(Ref(sois_dict), ks), "U" => "T"))
    sois.Motif = LongDNA{4}.(sois.Seq)
    sois.RCMotif = reverse_complement.(sois.Motif)
    sois.Index = 1:size(sois, 1)
    sois
end

"""
    count_mir_threads(mirseeds, seqs)

    Count exact matches of mirseeds in sequences (aimed at UTRs) using available threads.
    mirseeds and seqs are vectors of strings.
"""
function count_mir_threads(mirseeds, seqs)

    C = Matrix{Int}(undef, length(seqs), length(mirseeds))
    p = Progress(length(C), 0.5, "Counting mir seeds: ")
    Threads.@threads for ic in CartesianIndices(C)
        i, j = ic.I
        C[i, j] = count(mirseeds[j], seqs[i]) 
        next!(p)
    end
    finish!(p)
    C
end

"""
    collapse_isocounts(geneinds, iso_scores; collapse_fun=sum)

    Collapse isoform scores to gene scores using a given function (default: sum).

"""
function collapse_isoscores(geneinds, iso_scores::Matrix{T}; collapse_fun=sum) where {T}
   
    C = Matrix{T}(undef, length(geneinds), size(iso_scores, 2))
    
    p = Progress(length(geneinds), 0.5, "Collapsing iso counts: ")
    Threads.@threads for j = 1:size(iso_scores, 2)
        for i = 1:length(geneinds)
            C[i, j] = collapse_fun(iso_scores[geneinds[i], j])
        end
        next!(p)
    end
    finish!(p)
    C
end


"""
    match_utr_matrix(targetgenes, utrgenes, M::Matrix{T}; p = Progress(length(utrgenes))) where {T}

    Aligh a matrix M of counts in UTRs index by UTR genes to a list of target genes (targetgenes)
"""

function match_utr_matrix(targetgenes, utrgenes, M::Matrix{T}; p = Progress(length(utrgenes))) where {T}
    Z = zeros(T, length(targetgenes), size(M, 2)) 
    gi = Dict(tg => i for (i, tg) in enumerate(targetgenes))
    
    for (i, u) in enumerate(utrgenes)
        if haskey(gi, u)
            k = gi[u]
            Z[k, :] .= M[i, :]
        end
        next!(p)
    end
    Z
end

"""
    mir_count_dictdict
"""
function mir_count_dict(gex_pa_rz, utrs, gene_utrs, mir_iso_counts::Matrix{T}, mir_gene_counts::Matrix{T}) where {T}
   
    countdict = Dict{String, Matrix{T}}()
    p = Progress(length(gex_pa_rz.gpgenes), 0.1)
    Z_gene_tpm = match_utr_matrix(gex_pa_rz.gpgenes, gene_utrs.Gene, mir_gene_counts, p=p);
    finish!(p)
    
    countdict["pa"      ] = Z_gene_tpm
    countdict["rz"      ] = Z_gene_tpm
    

    countdict
end



"""
    count_seq_threads(motifs, seqs)

    Count sequences motifs in sequences using available threads.
    Motifs amd seqs are LongDNA{4} objects, allowing for ambiguity codes, with ExactSearchQuery
"""
function count_seq_threads(motifs, seqs)
    queries = ExactSearchQuery.(motifs, iscompatible)
    C = Matrix{Int}(undef, length(seqs), length(motifs))
    p = Progress(length(C), 0.5, "Counting motifs: ")
    Threads.@threads for ic in CartesianIndices(C)
        i, j = ic.I

        C[i, j] = length(findall(queries[j], seqs[i])) ## need to ensure that gaps are not respresented as N but as - for this to work
   
        next!(p)
    end
    finish!(p)
    C
end