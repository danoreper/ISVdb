source("./loadParamsFunc.R")
unlink(".RData")

###REPLACE_inFile###

inFile = parseCommandArgs()$inpArgs[1]
print(inFile)

x = try(load(inFile))
if(class(x)=="try-error")
{
    print(inFile)
    stop("failed to load input")
}

print(command)

## if(exists("prop"))
## {
##     for(nm in names(prop))
##     {
##         print(nm)
##         print(prop[[nm]])
##     }
## }


for(sourceFile in filesToSource)
{
    source(sourceFile)
}

load(funcFile)

## print("**********************************")
## print(prop$mnp$limit)
failed = c()
clusterOut = list()
for(i in 1:length(funcArgs))
{
    argz = c(funcArgs[[i]], otherGlobals)
    clusterOut[[i]] = try(do.call(func, argz))
    if(length(clusterOut[[i]])==1 && class(clusterOut[[i]])=="try-error")
    {
        ##print(clusterOut[[i]])
        failed = c(failed, i)
    }
}

x = try(save(list=c("clusterOut"), file=outFile))
y = try(save(list=c("failed"),     file=failFile))

if(length(failed)>0)
{
    print(command)
}
