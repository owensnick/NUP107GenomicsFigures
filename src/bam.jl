

function load_mir427_bam()
    bamfiles = glob("data/mir427/*.STAR.mir427.bam")
    bamfiles
end

function mir427_locus()
    chrom= "Chr03"
    loc = 133185000:133320000

    gaps = [133242393  133242665
            133245153  133250063
            133257403  133257678]

    chrom, loc, gaps
end

function mir427_locus_figure(gex_pa_rz, mir427_rz_pile)


    
    H_UC = mir427_rz_pile.H[:, gex_pa_rz.meta_rz.Treat .== "UC"]
    H_MO = mir427_rz_pile.H[:, gex_pa_rz.meta_rz.Treat .== "MO"]

    loc = mir427_rz_pile.loc
    dt = mir427_rz_pile.dt 
    pu = heatmap(first(loc):dt:last(loc), gex_pa_rz.meta_rz.Time[gex_pa_rz.meta_rz.Treat .== "UC"], log10.(H_UC' .+ 1), size=(1000, 100), clims=(0, 4), colorbar_ticks=0:4)
    pm = heatmap(first(loc):dt:last(loc), gex_pa_rz.meta_rz.Time[gex_pa_rz.meta_rz.Treat .== "UC"], log10.(H_MO' .+ 1), size=(1000, 100), clims=(0, 4), colorbar_ticks=0:4)
    
    plotscale = 1e+3
    t_uc = gex_pa_rz.meta_rz.Time[gex_pa_rz.meta_rz.Treat .== "UC"]
    y_uc = mir427_rz_pile.totalreads[gex_pa_rz.meta_rz.Treat .== "UC"]./plotscale
    
    t_mo = gex_pa_rz.meta_rz.Time[gex_pa_rz.meta_rz.Treat .== "MO"]
    y_mo = mir427_rz_pile.totalreads[gex_pa_rz.meta_rz.Treat .== "MO"]./plotscale
    st = range(first(t_uc), last(t_uc), length=100)
    gp_uc = gp_reg(t_uc, y_uc, st)
    gp_mo = gp_reg(t_mo, y_mo, st)
    
    
    pt = scatter(gex_pa_rz.meta_rz.Time, mir427_rz_pile.totalreads./plotscale, group=gex_pa_rz.meta_rz.Treat, grid=false)
    pt = plot(grid=false, ylabel="Normalised reads (millions)", xlabel="Time (hpf)")
    scatter!(t_uc, y_uc, c=:steelblue, lab="",marker=(stroke(0), 3),)
    plot!(st, gp_uc.sf, ribbon = [gp_uc.sf .- gp_uc.cil gp_uc.sf .- gp_uc.sf], lab="UC", c=:steelblue, fillalpha=0.2)
    
    # scatter!(t_mo, y_mo, c=:steelblue, lab="", markersize=2)
    scatter!(t_mo, y_mo, c=:orange, lab="",marker=(stroke(0), 3))

    plot!(st, gp_mo.sf, ribbon = [gp_mo.sf .- gp_mo.cil gp_mo.sf .- gp_mo.sf], lab="MO", c=:orange, fillalpha=0.2)

    
    plot!(pu, ylabel="UC", xformatter=x -> "", left_margin=5mm, bottom_margin=-2mm)
    plot!(pm, ylabel="MO", xformatter=x -> "", left_margin=5mm, bottom_margin=-2mm, top_margin=-1mm)
    
    
    xl_loc = (first(loc), last(loc))
    
    _, _, locgaps = mir427_locus()
    pg = plot(colorbar=true, ylims=(0.5, 1.5), grid=false, yticks=false, axis=false, top_margin=-1mm)
    for r in eachrow(locgaps)        
        plot!(rectangle(r[2] - r[1], 0.5, r[1], 0.75), zcolor=1, c=:black, lab="", line=stroke(0))
    end
    plot!(rectangle(100_000, 0.5, loc[1] + 20_000, 0.75), c=:black, lab="")
    
    plot!(pu, xlims=xl_loc)
    plot!(pm, xlims=xl_loc)
    plot!(pg, xlims=xl_loc)
        
    layout = @layout [[a ; b ; d{0.05h}] c{0.3w}]
    
    
    plot(pu, pm, pg, pt, layout=layout, size=(800, 200), fontfamily="helvetica")
end
rectangle(w, h, x, y) = Shape(x .+ [0,w,w,0], y .+ [0,0,h,h])

function mir427_pile(gex_pa_rz; recalc=false, dt=20, normscale=1e+6, mqs=-1, save=true, projdir=getprojectdir())
    bamfiles = joinpath.(dirname.(gex_pa_rz.meta_rz.File), "mir427", replace.(basename.(gex_pa_rz.meta_rz.File), ".genes.results" => ".STAR.mir427.bam"));
    chrom, loc, gaps = mir427_locus()
    pilefile = joinpath(projdir, "data", "mir427_rz_pile_$(chrom)_$(first(loc))_$(last(loc)).tsv.gz") 
    totalreadsfile = joinpath(projdir, "data", "mir427_rz_totalreads_$(chrom)_$(first(loc))_$(last(loc)).tsv.gz")
               
    

    if recalc && all(isfile, bamfiles)
        bps = @showprogress [ambigbampile(chrom, loc, bf, dt=dt, mqs=mqs) for bf in bamfiles];
        totalreads = normscale*last.(bps)./gex_pa_rz.meta_rz.TotalCounts
        H = normscale*mapreduce(first, hcat, bps)./gex_pa_rz.meta_rz.TotalCounts'
        if save
            CSV.write(pilefile, DataFrame(H, gex_pa_rz.meta_rz.Label), delim='\t', compress=true)
            CSV.write(totalreadsfile, DataFrame(C=totalreads), delim='\t', compress=true)
        end
        return (; chrom, loc, dt, H, totalreads)

    else
        if !isfile(pilefile) || !isfile(totalreadsfile)
            error("Precalculated files not found: $(pilefile) or $(totalreadsfile), please supply bam files")
        end
        H = Matrix(CSV.read(pilefile, DataFrame))
        totalreads = CSV.read(totalreadsfile, DataFrame).C
        return (; chrom, loc, dt, H, totalreads)
    end

end


"""
    ambigbampile(chrom, loc, bamfile; dt=10, mqs=255)

    Pileup reads at a given location in a bam file. Designed for use with Xt miR427 locus that is very repetitive, 
"""
function ambigbampile(chrom, loc, bamfile; dt=10, mqs=255)
    reader = open(BAM.Reader, bamfile, index=string(bamfile, ".bai"))
    n = 0
    xp = first(loc):dt:last(loc)

    cov = zeros(length(xp))
   
    # totalreads = DataStructures.counter(String)

    reads = Dict{String, Vector{BAM.Record}}()

    for record in eachoverlap(reader, chrom, loc)
        n += 1
        if BAM.mappingquality(record) < mqs
            continue
        end
        !haskey(reads, BAM.tempname(record)) && (reads[BAM.tempname(record)] = Vector{BAM.Record}())
        push!(reads[BAM.tempname(record)], record)
    end
    totalreads = length(reads)

    for (name, records) in reads

        mr = BAM.mappingquality.(records)
        for record in records
            aln = BAM.alignment(record)
            for i = 1:BAM.seqlength(record)
                k, op = seq2ref(aln, i)
                if ismatchop(op) && (k ∈ loc)
                    ind = cld(k - first(loc) + 1, dt)
    
                    cov[ind] += 1/length(records) 
                end
            end
        end
    end
    
    close(reader)
    cov, totalreads
end


function seq_frag_sizes_plot( ; file = "nup107_ribozero_seqfragsizes.tsv.gz", projdir=getprojectdir())

    filepath = joinpath(projdir, "data", file)

    fs = CSV.read(filepath, DataFrame)

    totalfrags = combine(groupby(fs, [:ReadLength]), :Count => sum => :Count)
    @with totalfrags bar(:ReadLength, :Count/1e+6, fill=0, fillalpha=0.2, leg=false, line=stroke(0))
    vline!([22], c=:black, ls=:dash, lab="")
    p = plot!(grid=false, xlabel="Sequenced Fragment Length", ylabel="Total Fragments (millions)", fontfamily="helvetica")
    ip = bar!(totalfrags.ReadLength, totalfrags.Count/1e+6, xlims=(14, 30), xticks=14:2:30, grid=false, ylims=(0, 0.1), fillalpha=0.2, leg=false, line=stroke(0), inset=(1, bbox(.085, .45, .35, .5, :bottom, :left)), subplot=2, fontfamily="helvetica")
    plot!(size=(400, 300))
end

"""

    load_short_seq_library(dir=""; recalc=true, file="nup107_ribozero_shortseqs.tsv.gz", projdir=getprojectdir())

    Loads file of short sequences (those shorter than the fragment size), collates from dir output if necessay saving those with size <= `n`
"""
function load_short_seq_library(dir=""; n=30, recalc=true, file="nup107_ribozero_shortseqs_le30.tsv.gz", projdir=getprojectdir())
    filepath = joinpath(projdir, "data", file)
    if isdir(dir) && recalc
        files = glob("*.gz", dir);
        dfs = CSV.read.(files, DataFrame)
        labels = first.(split.(basename.(files), "."));
        [df[!, :Label] .= l for (df, l) in zip(dfs, labels)]
        
        df = sort!(@subset(reduce(vcat, dfs), length.(:Seq) .<= n), :Seq)
        CSV.write(filepath, df, delim='\t', compress=true)
    else
        if !isfile(filepath)
            erorr("$filepath not found, rerun setting dir to regenerate")
        end
        df = CSV.read(filepath, DataFrame)
    end
    df.Length = length.(df.Seq)

    df
end

function annotate_mir427_counts!(df; mir427_seq="GAAAGUGCUUUCUGUUUUGGGCG", mm=4)
    mir427 = LongDNA{4}(replace(mir427_seq, "U" => "T"))
    mir427_rc = reverse_complement(mir427)
    df.FC = count.(Ref(string(mir427)), df.Seq) # exact match on forward strand
    df.RC = count.(Ref(string(mir427_rc)), df.Seq); # exact match on reverse
    fquery = ApproximateSearchQuery(mir427)
    rquery = ApproximateSearchQuery(mir427_rc)
    df.FA = occursin.(Ref(fquery), mm, LongDNA{4}.(df.Seq)); # approx match on forward up to mm mismatches
    df.RA = occursin.(Ref(rquery), mm, LongDNA{4}.(df.Seq)); # approx match on reverse up to mm mismatches

    df
end