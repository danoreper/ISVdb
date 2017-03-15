# 
# Loads params from yaml file.
###############################################################################
source("./loadParamsFunc.R")

prop = loadParams()
print("done loading")
for (aname in ls(prop))
{
    print("*********")
    print(aname)
    print(prop[[aname]])
}

##intentionally in global scope for ease of typing as its everywhere
dat = function(x)
{
    file.path(prop$data, x)
}

fp = file.path
