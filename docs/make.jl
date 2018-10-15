using Discuit
using Documenter

# create docs
makedocs(modules=[Discuit], doctest=true)
# deploy to GitHub
deploydocs(deps = Deps.pip("mkdocs", "python-markdown-math", "mkdocs-bootswatch") #"mdx_bib", "mkdocs-bootstrap386"
    , repo = "github.com/mjb3/Discuit.jl.git"
    , julia  = "1.0"
    , osname = "linux")
