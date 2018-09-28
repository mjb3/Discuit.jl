using Discuit
using Documenter
# create docs
makedocs(doctest = false)
# deploy github
deploydocs(repo = "github.com/mjb3/Discuit.jl.git")

## pdf (didn't work)
# using DocumenterLaTeX
# makedocs(format = :latex)
