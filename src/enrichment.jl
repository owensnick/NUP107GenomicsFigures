
## hypergeometric test/ Fisher's exact test
function ht(indA, indB)
    a = sum(indA .& indB)
    b = sum(indA .& .!indB)
    c = sum(.!indA .& indB)
    d = sum(.!indA .& .!indB)
    
    ft = FisherExactTest(a, b, c, d)
    or = (a/b)/(c/d)
    logor_se = sqrt(1/a + 1/b + 1/c + 1/d)
    cil = exp(log(or) - 1.96*logor_se)
    ciu = exp(log(or) + 1.96*logor_se)
    pvl = pvalue(ft, tail=:left)
    pvr = pvalue(ft, tail=:right)
    pvb = pvalue(ft, tail=:both)
    
    (; a, b, c, d, or, logor_se, cil, ciu, pvl, pvr, pvb)
    
end

function mirenrich(clusters, mirs, M, t = 1)
    uc = unique(clusters) 
    me = @showprogress [(Cluster=u, MirIndex=mirs[i], ht(clusters .== u, M[:, i] .>= t)...) for u in uc, i = 1:length(mirs)]
    DataFrame(vec(me))
end