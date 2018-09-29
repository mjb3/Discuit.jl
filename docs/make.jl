using Discuit
using Documenter

# create docs
makedocs(modules=[Discuit], doctest=true)
# deploy to GitHub
deploydocs(deps = Deps.pip("mkdocs", "python-markdown-math", "mkdocs-bootstrap386")
    , repo = "github.com/mjb3/Discuit.jl.git"
    , julia  = "nightly"
    , osname = "linux")

## pdf (didn't work)
# using DocumenterLaTeX
# makedocs(format = :latex)
