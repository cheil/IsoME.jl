using Documenter, DocumenterCitations, isoME   

DocMeta.setdocmeta!(
    isoME,
    :DocTestSetup,
    :(using isoME),
    recursive = true,
)

bib = CitationBibliography(
    joinpath(@__DIR__, "src", "refs.bib"),
    style = :authoryear,
)

makedocs(
    sitename = "IsoME.jl",
    modules = [isoME],
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing)=="true",
    ),
    pages = [
        "Home" => "home.md",
    ],
    warnonly = false,
    doctest = true,
    plugins = [bib],
)

deploydocs(
    repo = "https://gitlab.tugraz.at/computational-materials-group/codes/isotropic-me",
    push_preview = true,
)