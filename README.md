# Nup107 Genomics Figure Notebook


This julia repository contains data and code to perform analysis and generate figures for:

Valentyna Kostiuk, Rakib Kabir, Kimberly Morgan, Nick D. L. Owens, C. Patrick Lusk, Mustafa K. Khokha. *Nup107 contributes to the maternal to zygotic transition by preventing the premature nuclear export of pri-miRNA 427*

Repository contains the source data and code to understand the transcriptional response to depletion of *nup107* in early *Xenopus* development.

## Prerequistes
Julia >= 1.10, all julia packages and their versions are specified in the included Project.toml and Manifest.toml.

Additional this package uses https://github.com/exeter-tfs/MotifScanner.jl, which in turn employs Python matplotlib through https://github.com/JuliaPy/PyCall.jl.

## Installation
```bash
git clone https://github.com/owensnick/NUP107GenomicsFigures.jl
cd NUP107GenomicsFigures.jl
julia
```

Within julia activiate the current directory
```julia
] # to enter into Pkg mode
activate .
instantiate ## for first time installation
```

To regenerate figures either use jupyter notebook within `notebooks` directory or use script as follows:
```julia
 include("notebooks/nup107_figures.jl")
 ```
This will generate a `figs` folder and will generate all figure panels in `svg` format.