using Zubkov3
using Documenter

DocMeta.setdocmeta!(Zubkov3, :DocTestSetup, :(using Zubkov3); recursive=true)

makedocs(;
    modules=[Zubkov3],
    authors="Evgeny <e.a.malitsky@gmail.com'> and contributors",
    sitename="Zubkov3.jl",
    format=Documenter.HTML(;
        canonical="https://Vldjn.github.io/Zubkov3.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/Vldjn/Zubkov3.jl",
    devbranch="main",
)
