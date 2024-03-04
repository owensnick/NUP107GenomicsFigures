# setup project and load packages

using ThreadsX
using Markdown
using DataFrames, DataFramesMeta, CSV, Glob, ProgressMeter, CodecZlib
using Plots, StatsPlots, Measures
using Statistics, StatsBase, MultivariateStats, HypothesisTests
using LinearAlgebra, GaussianProcesses, Distributions, Optim, LineSearches, Suppressor
using Distances, Random, Clustering, ClusterOrderTools
using FASTX, BioSequences, Kmers, MotifScanner
using GenomicFeatures, BioAlignments, XAM

import DataStructures


include("rnaseq.jl")
include("utrs.jl")
include("plots.jl")
include("gp.jl")
include("clustering.jl")
include("enrichment.jl")
include("motifs.jl")
include("bam.jl")

function showwide(table)
    if !haskey(ENV, "DATAFRAMES_COLUMNS")
        ENV["DATAFRAMES_COLUMNS"] = "100"
    end
    c = ENV["DATAFRAMES_COLUMNS"]
    ENV["DATAFRAMES_COLUMNS"] = "10000"
    display(table)
    ENV["DATAFRAMES_COLUMNS"] = c;
    nothing;
end

function showwl(rows=200, cols=10000)
    function showwl(table)
        if !haskey(ENV, "DATAFRAMES_COLUMNS")
            ENV["DATAFRAMES_COLUMNS"] = "100"
        end
        if !haskey(ENV, "DATAFRAMES_ROWS")
            ENV["DATAFRAMES_ROWS"] = "25"
        end
        c = ENV["DATAFRAMES_COLUMNS"]
        r = ENV["DATAFRAMES_ROWS"]
        ENV["DATAFRAMES_COLUMNS"] = string(cols)
        ENV["DATAFRAMES_ROWS"] = string(rows)
        # display(table[1:min(rows, size(table, 1)), :])
        display(table)
        ENV["DATAFRAMES_COLUMNS"] = c;
        ENV["DATAFRAMES_ROWS"] = r;
        nothing;
    end
    
end

function getprojectdir()
    d = pwd()
    if basename(d) == "notebooks"
        return dirname(d)
    else
        return d
    end
end

function markdowntitle(x)
    md"""### $x
    ---""" |> display
end
