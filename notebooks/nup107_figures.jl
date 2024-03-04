ENV["GKSwstype"] = "100" ### this needed for GR when running on server without a display terminal, comment out should you want GR to display figures rather than saving

include(joinpath("..", "src", "project.jl"))

nups = loadnupgenes()
nupgeneplot(nups, "PA")
savedisplay("nup_gene_xep_pa.svg", dir="draft")
nupgeneplot(nups, "RZ")
savedisplay("nup_gene_xep_rz.svg", dir="draft");

gex_pa_rz = load_tpm_pa_rz();

writetpm(gex_pa_rz);

pc_tpm = pa_rz_cor_heatmap(gex_pa_rz, field=:tpm, dt=7,  title="TPM")
plot!(size=(480, 380))
savedisplay("corr_plot_parz.svg", dir="draft")

ppa = pcaplot(gex_pa_rz, parz="pa", xlims=(-50, 70))
plot!(size=(300, 200))
savedisplay("pca_pa.svg", dir="draft")
ppa = pcaplot(gex_pa_rz, parz="rz", dir=-1, xlims=(-50, 70))
plot!(size=(300, 200))
savedisplay("pca_rz.svg", dir="draft");

GP_PA = @suppress_out gp_genes_dt(gex_pa_rz.f_ind, gex_pa_rz.meta_pa, gex_pa_rz.tpm_pa; neldermead=true, threaded=true);
GP_RZ = @suppress_out gp_genes_dt(gex_pa_rz.f_ind, gex_pa_rz.meta_rz, gex_pa_rz.tpm_rz; neldermead=true, threaded=true);

gp_params_pa = paramtable(gex_pa_rz.gpgenes, GP_PA)
gp_params_rz = paramtable(gex_pa_rz.gpgenes, GP_RZ);


## build gp meta and tpm table for convenient plotting/clustering
pa_gp = gp_tpm_table(gex_pa_rz.tpm_pa, GP_PA, gex_pa_rz.filtind);
rz_gp = gp_tpm_table(gex_pa_rz.tpm_rz, GP_RZ, gex_pa_rz.filtind);
gpdata = Dict("pa" => (pa_gp.meta, pa_gp.tpm), "rz" => (rz_gp.meta, rz_gp.tpm)) ;

markdowntitle("PA/RZ tpm/gene\nDiagnostic plot of lengthscales in full and individual models")

p_ell_pa = @with gp_params_pa histogram2d(log.((:ℓ_u .+ :ℓ_m)./2), log.(:ℓ_f), bins=100, title="PolyA+")
p_ell_rz = @with gp_params_rz histogram2d(log.((:ℓ_u .+ :ℓ_m)./2), log.(:ℓ_f), bins=100, title="RZ")

map(p -> plot!(p, [0, 7], [0, 7], lab="", c=:red, xlabel="ℓ Mean U + M", ylabel="ℓf"), [p_ell_pa, p_ell_rz])
map(p -> plot!(p, [0, 7], [0, 7] .- 2, ls=:dash, lab="", c=:red), [p_ell_pa, p_ell_rz])
map(p -> plot!(p, [0, 7], [0, 7] .+ 2, ls=:dash, lab="", c=:red), [p_ell_pa, p_ell_rz])

    
plot(p_ell_pa, p_ell_rz, size=(800, 300), margin=5mm) |> display
plot(gp_param_hist(GP_PA), gp_param_hist(GP_RZ), size=(800, 300)) |> display

gp_de_all = parz_de_table(gex_pa_rz.gpgenes, GP_PA, GP_RZ);
markdowntitle("PA/RZ Gene TPM")
plot(plot_de_summary.(Ref(gp_de_all), ["PA", "RZ"])..., size=(1200, 250)) |> display
markdowntitle("PA/RZ Gene TPM"); 
lr_cd_histogram_parz(gp_de_all) |> display

parz_de_summary(gp_de_all)
plot!(size=(300, 300))
savedisplay("LR_parz.svg", dir="draft");

clinds = Dict{String, BitArray{1}}()
clinds["pa"]         = @with gp_de_all :PA_sig .& .!occursin.("ERCC-", :Gene);
clinds["rz"]         = @with gp_de_all :RZ_sig .& .!occursin.("ERCC-", :Gene);

ks = sort(collect(keys(clinds)))

cldf = DataFrame(K=ks, DataSet=first.(split.(ks, "_")), Condition=[join(split(k, "_")[2:end], "_") for k in ks], Count = [sum(clinds[k]) for k in ks])
unstack(cldf, :Condition, :DataSet, :Count)

theme(:wong)
markdowntitle("K means clustering silhouette scores for PA")
kmeans_ud_silhouette(pa_gp.meta, pa_gp.tpm, clinds["pa"], seed=1618033) |> display
markdowntitle("K means clustering silhouette scores for RZ")
kmeans_ud_silhouette(pa_gp.meta, pa_gp.tpm, clinds["rz"], seed=1618033) |> display

theme(:wong2)

# clustering function that allows for multiple ks and sets of DE indicioes
ks = [3]
kls = ["pa", "rz"]
p = Progress(length(kls)*length(ks))
KMKD = Dict(string(kl, "_", k) => (next!(p); kmeans_ud(gpdata[first(split(kl, "_"))]..., clinds[kl], k, k, seed=1618033)) for kl in kls, k in ks);
finish!(p);

parz_clusters = join_cluster_tables(KMKD["pa_3"], KMKD["rz_3"])
@assert gp_de_all.Gene == parz_clusters.Gene
parz_clusters = [parz_clusters gp_de_all[!, r"(PA|RZ)_(LR|CD)"]]
combine(groupby(parz_clusters, names(parz_clusters, r"PA_Cluster|PA_[UD]_ind")), nrow => :count) |> display

pa_up = @with parz_clusters sum((:PA_Ind .& :PA_U_ind))
pa_dn = @with parz_clusters sum((:PA_Ind .& .!:PA_U_ind))
rz_up = @with parz_clusters sum((:RZ_Ind .& :RZ_U_ind))
rz_dn = @with parz_clusters sum((:RZ_Ind .& .!:RZ_U_ind))
groupedbar(["polyA+", "ribozero"], [pa_up pa_dn ; rz_up rz_dn], lab=["Increased" "Decreased"], c=[:steelblue :orange], ylabel="Total Genes")

annotate!((0.5-.2, pa_up, text(pa_up, font(:bottom, "helvetica", 10))))
annotate!((0.5+.2, pa_dn, text(pa_dn, font(:bottom, "helvetica", 10))))
annotate!((1.5-.2, rz_up, text(rz_up, font(:bottom, "helvetica", 10))))
annotate!((1.5+.2, rz_dn, text(rz_dn, font(:bottom, "helvetica", 10))))
plot!(size=(200, 300), grid=false, fontfamily="helvetica", leg=:outertop)
savedisplay("total_up_down.svg", dir="draft");

markdowntitle("Total Agreements in differentially expressed genes U/D between PA and RZ")

parz_clusters.PAC = first.(parz_clusters.PA_Cluster)
parz_clusters.RZC = first.(parz_clusters.RZ_Cluster)

parz_da = unstack(combine(groupby(parz_clusters, [:PAC, :RZC]), nrow => :count), :PAC, :RZC, :count)
M = Matrix(parz_da[!, 2:end])
E = sum(M, dims=2)*sum(M, dims=1)/sum(M)
P = abs.(-log10.(ccdf.(Poisson.(E), M)))
F = M./E
M[1, 1] = 0
E[1, 1] = 0
P[1, 1] = 0
pm = heatannot(names(parz_da)[2:end], parz_da.PAC, M, size=(140, 140), colorbar=false, c=:blues, fs=font("helvetica", 9, :white))
savedisplay("intersection_parz_counts.svg", dir="draft")
pe = heatannot(names(parz_da)[2:end], parz_da.PAC, E, size=(140, 140), colorbar=false, c=:blues, fs=font("helvetica", 9, :white))
savedisplay("intersection_parz_expected.svg", dir="draft")
pp = heatannot(names(parz_da)[2:end], parz_da.PAC, P, size=(140, 140), colorbar=false, trf=x -> log10(x .+ 1),  c=:blues, fs=font("helvetica", 9, :white))
savedisplay("intersection_parz_poisson.svg", dir="draft");

markdowntitle("Clusters pa_3")

plotclusterud_stack(KMKD["pa_3"])
plot!(size=(600, 250), left_margin=5mm)
savedisplay("full_clusters_pa_3.svg", dir="draft")

markdowntitle("Clusters rz_3")
plotclusterud_stack(KMKD["rz_3"], lta = @layout [a ; b{0.15h}])
plot!(size=(600, 250), left_margin=5mm)
savedisplay("full_clusters_rz_3.svg", dir="draft");


parz_clstats = @chain parz_clusters begin
    groupby([:PA_Cluster, :RZ_Cluster])
    combine(nrow => :count)
    unstack(:RZ_Cluster, :PA_Cluster, :count)
    coalesce.(0)
    sort(:RZ_Cluster, by = clusterorder)
    
end
parz_clstats = parz_clstats[!, ["RZ_Cluster" ; sort(names(parz_clstats)[2:end], by=clusterorder)]]
M = Matrix(parz_clstats[!, 2:end])
E = sum(M, dims=2)*sum(M, dims=1)/sum(M)
P = abs.(-log10.(ccdf.(Poisson.(E), M)))
F = M./E
M[1, 1] = 0
E[1, 1] = 0
F[1, 1] = 0
pc = heatannot(names(parz_clstats)[2:end], parz_clstats.RZ_Cluster, M, c=:blues, fs=font("helvetica", 8, :white)); annotate_parz_clusterplot!()
plot!(title="Intersection between PA and RZ clusters")
pp = heatannot(names(parz_clstats)[2:end], parz_clstats.RZ_Cluster, P, c=:blues, fs=font("helvetica", 8, :white)); annotate_parz_clusterplot!()
plot!(title="-log10 Poisson P-value\noverrepresentation")
plot(pc, pp, size=(610, 290), margin=5mm, fontfamily="helvetica", titlefont=font(12, "helvetica"), colorbar=false)
savedisplay("parz_cluster_intersect_v2.svg", dir="draft");


phu = viscluster_double(KMKD["pa_3"].meta.Time[KMKD["pa_3"].meta.Treat .== "UC"], KMKD["pa_3"].KMK_U, c=:viridis, colorbar=true, ylabel="Up", titlefont=font(10, "helvetica"), xlabel="Time", bottom_margin=3mm, title="UC  MO")
plot!(phu, size=(200, 375), fontfamily="helvetica") 
savedisplay("up_cluster_heatmap.svg", dir="draft")

p, phs = plotclustermeansd_vert(KMKD["pa_3"].KMK_U, KMKD["pa_3"].meta)    
plot!(phs[1], title="Cluster Mean", fontfamily="helvetica")
# plot!(phs[1], bottom_margin=-1mm, xticks=false)
[plot!(p, top_margin=-1mm, bottom_margin=-1mm,  xticks=false) for p in phs[1:end-1]]
plot!(phs[end], top_margin=-1mm, xlabel="Time")
plot!(phs[end],xlabel="Time", xticks=2.5:2.5:10)

pcs = plot(phs..., layout=(3, 1), size=(150, 350))
savedisplay("cluster_means.svg", dir="draft")


exgenes = ["hivep1", "supt16", "wisp3"]

ghs = [plotgpset(g, gex_pa_rz.gpgenes, GP_PA, showparams=false, samples=["U", "M"], fontfamily="helvetica", marker=(2, stroke(0))) for g in exgenes]
[plot!(p, top_margin=-1mm, bottom_margin=-1mm,  xticks=false) for p in ghs[1:end-1]]
plot!(ghs[end], top_margin=-1mm, xlabel="Time")
plot!(ghs[end],xlabel="Time")
[plot!(p, ylabel=g, title="") for (p, g) in zip(ghs, exgenes)]
plot!(ghs[1], title="Example Genes")
pgs = plot(ghs..., layout=(3, 1), size=(150, 350), framestyle=:box, grid=false)
savedisplay("example_genes.svg", dir="draft");


mirseqs, mirs   = loadmirs(organism=r"tropicalis");
utrs, gene_utrs = loadutrs(gex_tested=[gex_pa_rz]);
mir_iso_counts  = count_mir_threads(mirs.RCSeed, utrs.seq);
mir_gene_counts = collapse_isoscores(gene_utrs.Index, mir_iso_counts, collapse_fun=sum); # counts of mir seeds in utrs against utrs present in anntation file
mir_counts = mir_count_dict(gex_pa_rz, utrs, gene_utrs, mir_iso_counts, mir_gene_counts) # counts of mir seeds in utrs of genes quantified

miren_pa = @subset(innerjoin(mirs, mirenrich(KMKD["pa_3"].gtable.Cluster, mirs.MirIndex, mir_counts["pa"], 1), on=:MirIndex), :Cluster .!= "N")
miren_rz = @subset(innerjoin(mirs, mirenrich(KMKD["rz_3"].gtable.Cluster, mirs.MirIndex, mir_counts["rz"], 1), on=:MirIndex), :Cluster .!= "N");

### Motifs of interest
mois = load_motifs_of_interest()
mois_iso_counts = count_seq_threads(mois.Motif, LongDNA{4}.(utrs.seq));
mois_gene_counts = @time collapse_isoscores(gene_utrs.Index, mois_iso_counts, collapse_fun=sum);
mois_counts = mir_count_dict(gex_pa_rz, utrs, gene_utrs, mois_iso_counts, mois_gene_counts)

moten_pa = @subset(innerjoin(mois, mirenrich(KMKD["pa_3"].gtable.Cluster, mois.Index, mois_counts["pa"], 1), on=:Index => :MirIndex), :Cluster .!= "N")
moten_rz = @subset(innerjoin(mois, mirenrich(KMKD["rz_3"].gtable.Cluster, mois.Index, mois_counts["rz"], 1), on=:Index => :MirIndex), :Cluster .!= "N");

markdowntitle("3' UTR Motif enrichment")
plot(utr_motif_heatmap(moten_pa, clims=(0, 6), title="PA+ Cluster Enrichment"),
     utr_motif_heatmap(moten_rz, clims=(0, 6), title="RZ Cluster Enrichment"), layout=(1, 2), size=(750, 150), colorbar=false, fontfamily="helvetica")
savedisplay("utr_motif_enrichment.svg", dir="draft");

p_pa = mir_enrich_heatmap(miren_pa, title="PA+ Cluster Enrichment", fs=font("helvetica", 6))
p_rz = mir_enrich_heatmap(miren_rz, title="RZ Cluster Enrichment", fs=font("helvetica", 6))

plot(p_pa, p_rz, size=(900, 450), fontfamily="helvetica", colorbar=false, clims=(0, 7.5), ytickfont=font(7, "helvetica"), titlefont=font("helvetica", 12))
savedisplay("supp_mir_enrichment.svg", dir="draft");

markdowntitle("Do mir427 and ARE motifs occur in the same UTRs?")
mir427_ind = findfirst(mirs.RCSeed .== "GCACTT")
are_ind = findfirst(mois.ID .== "ARE")
selind = KMKD["pa_3"].gtable.Cluster .∈ Ref(["U1"])
# selind .= true
ht(mir_counts["pa"][selind, mir427_ind] .> 0, mois_counts["pa"][selind, are_ind] .> 0)

markdowntitle("Do mir427 and YTHDF2 motifs occur in the same UTRs?")
mir427_ind = findfirst(mirs.RCSeed .== "GCACTT")
yth_ind = findfirst(mois.ID .== "YTHDF2")
selind = KMKD["pa_3"].gtable.Cluster .∈ Ref(["U1"])
# selind .= true
ht(mir_counts["pa"][selind, mir427_ind] .> 0, mois_counts["pa"][selind, yth_ind] .> 0)

markdowntitle("Do ARE and YTHDF2 motifs occur in the same UTRs?")
selind = KMKD["pa_3"].gtable.Cluster .∈ Ref(["U1", "U2"])
# selind .= true
ht(mois_counts["pa"][selind, are_ind] .> 0, mois_counts["pa"][selind, yth_ind] .> 0)

proms, gene_proms = loadpromoters();
motsel = loadmotifselection();

mot_sel_enrich = mot_score_thresh_enrich(motsel, proms, gene_proms, parz_clusters, minscore=5.0);

phs = Plots.Plot[]
for gdf in groupby(mot_sel_enrich, [:MotifName, :MotifID])
   # display(gdf)
    
    labels = join.(gdf.ClusterGroup, "\n")
    labels[end] = "All"
    p = bar(1:size(gdf, 1), -log10.(gdf.pvr), c=1:size(gdf, 1), xticks=(1:size(gdf, 1), labels), title=string.(gdf.MotifName[1]), lab="")
    hline!([-log10(0.05)], c=:red, ls=:dash, lab="")
    plot!(right_margin=-5mm, ylims=(0, 5), xtickfont=font(6, "helvetica"))
    push!(phs, p)
    
end

ghs = [plotgpset(Regex(string("\\|", lowercase(m), "\$")), gex_pa_rz.gpgenes, GP_PA) for m in replace.(getindex.(motsel, :name), "Arid3a" => "arid4b")]
[plot!(p, title=m) for (p,m) in zip(ghs, replace.(getindex.(motsel, :name), "Arid3a" => "arid4b"))]
plot!(ghs[1], ylabel="TPM")
plot!(phs[1], ylabel="-log10 p")
plot(ghs..., phs..., layout=(2, length(phs)), left_margin=-2mm, right_margin=-2mm, size=(1100*0.75, 360*0.75),  fontfamily="helvetica", titlefont=font(8, "helvetica"), bottom_margin=5mm, grid=false)
savedisplay("motif_enrichment.svg", dir="draft");

mir427_rz_pile = mir427_pile(gex_pa_rz, mqs=0, recalc=false);


mir427_locus_figure(gex_pa_rz, mir427_rz_pile)
savedisplay("mir427_locus_figure.svg", dir="draft");

seq_frag_sizes_plot()
savedisplay("short_fragment_sizes.svg", dir="draft");

shortseqs = load_short_seq_library();
annotate_mir427_counts!(shortseqs, mm=4)
@with combine(groupby(shortseqs, :Length), :FC => sum => :FC, :RC => sum => :RC, :FA => sum => :FA, :RA => sum => :RA) bar(:Length, :RA)
xticks!(0:2:30)
yticks!(0:2:30)
plot!(xlabel="Fragment size", ylabel="Count", title="Total fragments containing mature mir427 microRNA sequence\nup to 4 mismatches", leg=false, fontfamily="helvetica", titlefont=font(12, "helvetica"))
