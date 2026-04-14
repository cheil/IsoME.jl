using Documenter, DocumenterCitations, IsoME   

DocMeta.setdocmeta!(
    IsoME,
    :DocTestSetup,
    :(using IsoME),
    recursive = true,
)

bib = CitationBibliography(
    joinpath(@__DIR__, "src", "refs.bib"),
    style = :authoryear,
)

makedocs(
    sitename = "IsoME.jl",
    modules = [IsoME],
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing)=="true",
    ),
    pages = [
        "Home" => "index.md",
        "Input" => "Input.md",
        "Best Practices" => "bestPractices.md",
        "FAQ"   => "FAQ.md",
    ],
    warnonly = false,
    doctest = true,
    plugins = [bib],
    checkdocs=:exports,
)

deploydocs(
    repo = "github.com/cheil/IsoME.jl",
    push_preview = true,
    versions = nothing,
    branch = "gh-pages",
    devbranch = "main",
)