using MagnetoIonic
using Documenter

DocMeta.setdocmeta!(MagnetoIonic, :DocTestSetup, :(using MagnetoIonic); recursive=true)

makedocs(;
    modules=[MagnetoIonic],
    authors="James Cannon",
    sitename="MagnetoIonic.jl",
    format=Documenter.HTML(;
        canonical="https://JamesMCannon.github.io/MagnetoIonic.jl",
        edit_link="master",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/JamesMCannon/MagnetoIonic.jl",
    devbranch="master",
)
