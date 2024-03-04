
dtm(m, α=1, β=1000) = x -> sqrt(β*(x/m) + α)
itm(m, α=1, β=1000) = y -> ifelse(y > sqrt(α), m*(y*y - α)/β, 0)
sqd(x) = x*x
function gp_reg(t, y, st;  datatrans=sqrt, invtrans=sqd, dty = datatrans.(y), initparams = (log(mean(dty)), 2.0, 2*log(std(dty))), f_prior = Normal(3, 2.0), ℓ_prior = Normal(2, 4.0), n_prior=Normal(0.0, 4), use_priors=true, neldermead=true)
    σf, ℓ, σn = initparams


    kern = Mat52Iso(ℓ, σf)

    use_priors && set_priors!(kern, [ℓ_prior, f_prior])
    gp = GP(t, dty, MeanZero(), kern, σn)
    
    use_priors && set_priors!(gp.logNoise, [n_prior] )

    try
        optimize!(gp) #method=LBFGS()) #HagerZhang is the default, sometimes throws an error switch to backtracking in that case
   catch err
        optimize!(gp, method=LBFGS(linesearch=LineSearches.BackTracking()))
    end
    neldermead && optimize!(gp, method=NelderMead())

  
    gp_pred(t, y, st, gp, invtrans=invtrans)
end

function gp_pred(t, y, st, gp; invtrans=sqd)
    μ, v = predict_y(gp, st)
    sf  = invtrans.(μ)
    cil = invtrans.(μ .- 1.96*sqrt.(v))
    ciu = invtrans.(μ .+ 1.96*sqrt.(v))
    (gp=gp, μ=μ, v=v, sf=sf, cil=cil, ciu=ciu, t=t, y=y, st=st)
end

"""
    gp_set_dt_time(t_u, t_m, stf, st, u_ind, yf, si; α=1, β=1000, neldermead=false)

    Set up a GP model for the UC and MO time series and predict the time series for the full set of time points.
"""
function gp_set_dt_time(t_u, t_m, stf, st, u_ind, yf, si; α=1, β=1000, neldermead=false)

    yu = yf[u_ind];
    ym = yf[.!u_ind];
    m = maximum(yf)

    dt = dtm(m, α, β)
    it = itm(m, α, β)


    gpru = gp_reg(t_u, yu, st, datatrans=dt, invtrans=it, neldermead=neldermead);
    gprm = gp_reg(t_m, ym, st, datatrans=dt, invtrans=it, neldermead=neldermead);

    σf = (gpru.gp.kernel.σ2 + gprm.gp.kernel.σ2)/2
    ℓ  = (gpru.gp.kernel.ℓ + gprm.gp.kernel.ℓ)/2
    σn = (exp(2*gpru.gp.logNoise.value) + exp(2*gprm.gp.logNoise.value))/2
    initparams = (log(σf), log(ℓ), log(σn))
        
    gprf = gp_reg(stf, yf[si], st, initparams=initparams, datatrans=dt, invtrans=it, neldermead=neldermead);


    #### test to see if the length scales for the full model are sufficiently different from individual models
    ul = gpru.gp.kernel.ℓ
    ml = gprm.gp.kernel.ℓ
    fl = gprf.gp.kernel.ℓ
    reset_f = false
    reset_u = false
    reset_m = false

    if fl > (ul + ml) # fl is twice the mean of ml and ul
        gprf.gp.kernel.ℓ = min(ul, ml)
        optimize!(gprf.gp)
        gprf = gp_pred(stf, yf[si], st, gprf.gp, invtrans=it)
        reset_f = true
        fl = gprf.gp.kernel.ℓ
    end

    if ul > (fl + ml) # ul is twice the mean of ml and fl
        gpru.gp.kernel.ℓ = min(fl, ml)
        optimize!(gpru.gp)
        gpru = gp_pred(t_u, yu, st, gpru.gp, invtrans=it)
        reset_u = true
        ul = gpru.gp.kernel.ℓ
    end

    if ml > (fl + ul) # ml is twice the mean of ul and fl
        gprm.gp.kernel.ℓ = min(fl, ul)
        optimize!(gprm.gp)
        gprm = gp_pred(t_m, ym, st, gprm.gp, invtrans=it)
        reset_m = true
        ml = gprm.gp.kernel.ℓ
    end
        
    (u=gpru, m=gprm, f=gprf, reset_u=reset_u, reset_m=reset=m, reset_f=reset_f)

end


"""
    gp_genes_dt(inds, meta, tpm; α=1, β=1000, neldermead=false, threaded=false, mapfun=ifelse(threaded, ThreadsX.map, map))

    Setup gp models for all genes and fit GP models potentially threaded.

"""

function gp_genes_dt(inds, meta, tpm; α=1, β=1000, neldermead=false, threaded=false, mapfun=ifelse(threaded, ThreadsX.map, map))
    

    u_ind = meta.Treat .== "UC"
    t_u = meta.Time[u_ind]
    t_m = meta.Time[.!u_ind]

    
    st = range(minimum(t_u), maximum(t_u), length=100)

    tf = [t_u ; t_m]
    si = sortperm(tf)
    stf = tf[si]

    T = Matrix(tpm[inds, meta.Label])'

    p = Progress(length(inds), 1, "GP: ")
    bt = BLAS.get_num_threads()
    threaded && BLAS.set_num_threads(1)
    GP = mapfun(i -> ( next!(p); gp_set_dt_time(t_u, t_m, stf, st, u_ind, T[:, i], si; α=1, β=1000, neldermead=neldermead)), 1:length(inds))
    threaded && BLAS.set_num_threads(bt)
    
    GP
end


### functions to calculate stats on GP regression

"""
    gp_bic(gprf, gprs)
    
    Calculate the Bayesian Information Criterion (BIC) for a full model and individual models for UC and MO.
"""
function gp_bic(gprf, gprs)
    f_mll = gprf.gp.mll
    s_mll = mapreduce(gpr -> gpr.gp.mll, +, gprs)
    lr = s_mll - f_mll
    bic = lr - 3*log(gprf.gp.nobs)/2
    lr, bic
end

"""
    maxcohensd(gpr, fieldA=:u, fieldB=:m)

    Calculate the maximum Cohen's d  between UC and MO.
"""
function maxcohensd(gpr, fieldA=:u, fieldB=:m)
    md = 0.0
    mdi = 0
    
    for (μA, μB, vA, vB) in zip(gpr[fieldA].μ, gpr[fieldB].μ, gpr[fieldA].v, gpr[fieldB].v)
        v = abs(μA - μB)/sqrt(0.5*(vA + vB))
        # v = abs(μA - μB)/sqrt(1.0*(vA + vB))
        if v > md
            md = v
            if μA > μB
                mdi = -1
            elseif μA < μB
                mdi = 1
            else
                mdi = 0
            end
        end
    end
    md*mdi
end

"""
    paramtable(genes, G, fields=[:u, :m, :f])

    Create a DataFrame with the parameters of the GP models for UC, MO and full model.
"""
function paramtable(genes, G, fields=[:u, :m, :f])
    σf = [[sqrt.(g[f].gp.kernel.σ2)      for g in G] for f in fields]
    ℓ  = [[g[f].gp.kernel.ℓ       for g in G] for f in fields]
    σn = [[exp.(g[f].gp.logNoise.value) for g in G] for f in fields]
    SNR = [f./n for (f, n) in zip(σf, σn)]

    mean_u = [mean(g.u.gp.y) for g in G]
    mean_m = [mean(g.m.gp.y) for g in G]
    mean_f = [mean(g.f.gp.y) for g in G]

    std_u = [std(g.u.gp.y) for g in G]
    std_m = [std(g.m.gp.y) for g in G]
    std_f = [std(g.f.gp.y) for g in G]
    
    


    df = DataFrame(σf, Symbol.(string.(:σf, "_", fields)))
    dl = DataFrame(ℓ, Symbol.(string.(:ℓ, "_", fields)))
    dn = DataFrame(σn, Symbol.(string.(:σn, "_", fields)))
    ds = DataFrame(SNR, Symbol.(string.(:SNR, "_", fields)))
    [DataFrame(Gene=genes, Index=1:length(genes)) df dl dn ds DataFrame(mean_u=mean_u, mean_m=mean_m, mean_f=mean_f) DataFrame(std_u=std_u, std_m=std_m, std_f=std_f)] 
end


"""
    gp_de_table(genes, GPC)

    Create a DataFrame with the log likelihood ratio, BIC and maximum Cohen's d for all genes.
"""
function gp_de_table(genes, GPC)
    lr_bic = [gp_bic(gp.f, (gp.u, gp.m)) for gp in GPC]
    DataFrame(Gene=genes, Index=1:length(genes), LR=first.(lr_bic), BIC=last.(lr_bic), CD=maxcohensd.(GPC))
end

"""
    parz_de_table(gpgenes, GP_PA, GP_RZ)

    Create dataframe of GP differential expression results for PA and RZ.
"""
function parz_de_table(gpgenes, GP_PA, GP_RZ)
    gpde_pa = gp_de_table(gpgenes, GP_PA)
    gpde_rz = gp_de_table(gpgenes, GP_RZ)
    @assert gpde_pa.Gene == gpde_rz.Gene

    gpde_all = [rename(gpde_pa, :LR => :PA_LR, :BIC => :PA_BIC, :CD => :PA_CD) rename(gpde_rz[!, [:LR, :BIC, :CD]], :LR => :RZ_LR, :BIC => :RZ_BIC, :CD => :RZ_CD)];
    gpde_all.PA_sig = gpde_all.PA_LR .> 0
    gpde_all.RZ_sig = gpde_all.RZ_LR .> 0;

    gpde_all.GeneName = last.(split.(gpde_all.Gene, "|", keepempty=false))
    gpde_all
end


function gp_tpm_table(tpm, GPU, exfiltind)
    dt = 3
    st = GPU[1].u.st[1:dt:end]
    SF = mapreduce(g -> mapreduce(f -> Float64.(g[f].sf[1:dt:end]), vcat, [:u, :m]), hcat, GPU);
    metagp = DataFrame(Treat=repeat(["UC", "MO"], inner=length(st)), Time=repeat(st, outer=2))
    metagp.Label = string.(metagp.Treat, "_", repeat(1:length(st), outer=2))
    tpm_gp = [DataFrame(Gene=tpm.Gene[exfiltind]) DataFrame(SF', Symbol.(metagp.Label))];
    
    (meta=metagp, tpm=tpm_gp)
end

##### Gaussian Process plotting functions



function gp_param_hist(G, fields=[:u, :m, :f]; f_prior= Normal(1.0, 3.0), ℓ_prior = Normal(1, 3.0), n_prior=Normal(0.5, 3.0))

    σf2    = [[g[f].gp.kernel.σ2      for g in G] for f in fields]
    ℓ      = [[g[f].gp.kernel.ℓ       for g in G] for f in fields]
    log_σn = [[g[f].gp.logNoise.value for g in G] for f in fields]

    
    σf_bins = range(0, 5, length=200)
    ℓ_bins  = range(0, 5, length=200)
    σn_bins = range(-3, 3, length=200)
    f_phs = [stephist(log.(sf)/2, bins=σf_bins, lab=string(l), fill=0, normed=true, c=:black, fillalpha=0.2, ylims=(0, 8), yticks=0:2:4) for (sf, l) in zip(σf2, fields)]
    ℓ_phs = [stephist(log.(el),   bins=ℓ_bins,  lab=string(l), fill=0, normed=true, c=:red, fillalpha=0.2, ylims=(0, 3), yticks=0:2) for (el, l) in zip(ℓ, fields)]
    n_phs = [stephist(sn,   bins=σn_bins, lab=string(l), fill=0, normed=true, c=:steelblue, fillalpha=0.2, ylims=(0, 2.5), yticks=0:2) for (sn, l) in zip(log_σn, fields)]

    # σf_bins = range(.5, 2, length=200)
    # f_phs = [stephist(log10.(sf)/2, bins=σf_bins, lab=string(l), fill=0, normed=true) for sf in σf2]

    ### get prios from f gp
    ℓ_prior = G[1].f.gp.kernel.priors[1]
    f_prior = G[1].f.gp.kernel.priors[2]
    n_prior = G[1].f.gp.logNoise.priors[1]

    for (i, (f, l, n)) in enumerate(zip(f_phs, ℓ_phs, n_phs))
        if i == 1
            plot!(f, title="sigma_f")
            plot!(l, title="ell")
            plot!(n, title="sigma_n")
        end
        if i < length(fields)
            plot!(f, xformatter=x -> "")
            plot!(l, xformatter=x -> "")
            plot!(n, xformatter=x -> "")
        end
        ylabel!(f, string(fields[i]))
        plot!(f, σf_bins, pdf.(f_prior, σf_bins))
        plot!(l, ℓ_bins, pdf.(ℓ_prior, ℓ_bins))
        plot!(n, σn_bins, pdf.(n_prior, σn_bins))
    end

    lt = @layout [[a ; b ; c] [d ; e ; f ] [g ; h ; i]]
    plot(f_phs..., ℓ_phs..., n_phs..., layout=lt, margin=-1mm, leg=false)

end


function plot_de_summary(detable, pattern="PA")
    
    bins = range(-20, 20, length=90)
    fields = ["Gene" ; names(detable, Regex(pattern))]
    
    
    gpde = @subset(rename!(detable[!, fields], replace.(fields, string(pattern, "_") => "")), .!occursin.(r"ERCC-", :Gene))
    # clind = gpde.LR .> 0
    ph = plot()
    stephist!(gpde.LR, fill=0, normed=false, fillalpha=0.2, bins=bins, c=[:steelblue :orange], lab=string("LR > 0: ", sum(gpde.sig)))
    vline!([0], c=:black, lab="", ls=:dash)
    plot!(xlabel="LR", title="$(pattern): Histogram of LR", ylabel="# Genes")
    
    pe = plot(leg=:topright, xlabel="LR", title="$(pattern): Inverse cumulative LR")
    plot!(bins, size(gpde, 1)*(1 .- ecdf(gpde.LR)(bins))/1000, c=:steelblue, lab=string("LR > 0: ", sum(gpde.sig)))
    vline!([0], c=:black, lab="", ls=:dash)
    hline!([sum(gpde.sig)/1000], c=:red, lab="")
    plot!(ylabel="# Genes (k)")
    annotate!([(20, (sum(gpde.sig))/1000, text(string(sum(gpde.LR .> 0), " genes"), font(:red, :right, :bottom, "helvetica", 9)))])
    plot(ph, pe, size=(800*1.25, 300*1.25), bottom_margin=5mm, left_margin=5mm, grid=false, fontfamily="helvetica", titlefont=font("helvetica", 10))
end

function parz_de_summary(detable; fix_xl=(0, 0.135))
    gpde = @subset(detable, .!occursin.(r"ERCC-", :Gene))
    total_de_summary = combine(groupby(gpde, [:PA_sig, :RZ_sig]), nrow => :count)
    
    @with gpde histogram2d(:PA_LR, :RZ_LR, bins=250)

    p = @with gpde marginalhist(:PA_LR, :RZ_LR, bins=250, density=true)
    vline!(p.subplots[1], [0],  c=:black, lab="", ls=:dash)
    vline!(p.subplots[2], [0],  c=:black, lab="", ls=:dash)
    hline!(p.subplots[2], [0],  c=:black, lab="", ls=:dash)
    hline!(p.subplots[3], [0],  c=:black, lab="", ls=:dash)
    plot!(p.subplots[3], xlims=fix_xl) ### manually set xlims on 3rd subplot as this is broken


    plot!(p.subplots[2], grid=false, xlabel="LR polyA+", ylabel="LR Ribozero")
    
    
    xl = xlims(p.subplots[2])
    yl = ylims(p.subplots[2])
    xd = xl[2] - xl[1]
    yd = yl[2] - yl[1]
    
    a, b, c, d = total_de_summary.count
    pv = pvalue(FisherExactTest(a, b, c, d), tail=:right)
    or = (a/b)/(c/d)

    c = [(row.PA_sig*.9*xd + xl[1] + .01*xd, row.RZ_sig*.9*yd + yl[1] + .05*yd, text(row.count, font("helvetica", 8, :left))) for row in eachrow(total_de_summary)]
    annotate!(p.subplots[2], c)
    plot!(p.subplots[2], [-25, 45], [-25, 45], c=:red, ls=:dash, lab="", title=string("p ", StatsBase.PValue(pv), " OR: ", tt(or)), titlefont=font("helvetica", 12))
    plot!(p.subplots[2], xlims=xl, ylims=yl)
end


function lr_cd_histogram_parz(gpde)
    p_pa_cd = @with gpde histogram2d(:PA_LR, :PA_CD, bins=200, title="PA")
    p_rz_cd = @with gpde histogram2d(:RZ_LR, :RZ_CD, bins=200, title="RZ")

    phs = [p_pa_cd, p_rz_cd]
    map(p -> hline!(p, [0, -2, 2], c=:black, lab="", ls=:dash), phs)
    map(p -> vline!(p, [0], c=:black, lab="", ls=:dash), phs)

    plot(phs..., link=:y, ylims=(-10, 10), size=(800, 300), xlabel="LR", ylabel="CD")
end
