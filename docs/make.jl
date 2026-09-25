using Documenter, IsoME

DocMeta.setdocmeta!(
    IsoME,
    :DocTestSetup,
    :(using IsoME),
    recursive = true,
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
        "Matsubara Solver" => "MatsubaraSolver.md",
        "Real Axis Solver" => "RealAxisSolver.md",
        "Best Practices" => "bestPractices.md",
        "Troubleshooting" => "Troubleshooting.md",
        "FAQ"   => "FAQ.md",
    ],
    warnonly = false,
    doctest = true,
    checkdocs=:exports,
)

deploydocs(
    repo = "github.com/cheil/IsoME.jl",
    push_preview = true,
    versions = nothing,
    branch = "gh-pages",
    devbranch = "main",
)