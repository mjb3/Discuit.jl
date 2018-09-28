using Discuit
using Documenter
# create docs
makedocs(doctest = false)
# deploy github
deploydocs(
    repo = "github.com/USER_NAME/PACKAGE_NAME.jl.git",
)

## pdf (didn't work)
# using DocumenterLaTeX
# makedocs(format = :latex)
