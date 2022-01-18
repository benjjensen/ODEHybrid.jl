using Documenter
using ODEHybrid

makedocs(
    sitename = "ODEHybrid.jl",
    format = Documenter.HTML(prettyurls = false),
    pages = [
        "Introduction" => "index.md",
        "API" => "api.md"
    ]
)

deploydocs(
    repo = "github.com/benjjensen/ODEHybrid.jl.git",
    devbranch = "main"
)