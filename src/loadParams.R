# 
# Script for loading params from yaml file.
###############################################################################
source("./loadParamsFunc.R")

prop = loadParams()
for (aname in ls(prop))
{
    ## print("*********")
    ## print(aname)
    ## print(prop[[aname]])
}

##intentionally in global scope for ease of typing as its everywhere
dat = function(...)
{
    file.path(prop$data, ...)
}

outf = function(...)
{
    apath = file.path(prop$output, ...)
    dir.create(dirname(apath), recursive = T, showWarnings = F)
    return(apath)
}
fp = file.path
rm("aname")
