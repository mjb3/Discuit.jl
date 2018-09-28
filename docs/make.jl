using Discuit
using Documenter
# create docs
# makedocs()
# deploy github
# deploydocs(repo = "github.com/mjb3/Discuit.jl.git", target = "docs/build/")

# create docs
makedocs(modules=[Discuit], doctest=true)

deploydocs(deps = Deps.pip("mkdocs", "python-markdown-math")
    , repo = "github.com/mjb3/Discuit.jl.git"
    , julia  = nightly
    , osname = "linux")

## pdf (didn't work)
# using DocumenterLaTeX
# makedocs(format = :latex)
