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

for(sourceFile in filesToSource)
{
    source(sourceFile)
}

load(funcFile)


failed = c()
clusterOut = list()
for(i in 1:length(funcArgs))
{
    clusterOut[[i]] = try(do.call(func, c(funcArgs[[i]], otherGlobals)))
    if(class(clusterOut[[i]])=="try-error")
    {
        failed = c(failed, i)
    }
}

x = try(save(list=c("clusterOut"), file=outFile))
y = try(save(list=c("failed"),     file=failFile))

