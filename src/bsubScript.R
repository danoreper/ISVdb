## A general purpose script that is used by any invocation of bsub paraallelization across the cluster.
##Loads in serialized function and function arguments, as well a file location to write to, calls the function, writes results to file
## func, funcArgs, and outFile must be defined in the session
source("./loadParamsFunc.R")

##prop = loadParams()

###REPLACE_inFile###

inFile = parseCommandArgs()$inpArgs[1]
print(inFile)

counter = 1
##sometimes the file system fails due to large numbers of temp files being created... keep trying to reload the file a few times with a pause in between attempts.
while(counter<20)
{
    x = try(load(inFile))
    if(class(x)=="try-error")
    {
        print(paste0("trying counter ", counter))
        Sys.sleep(min(abs(50+rnorm(0,5)), 120))
        counter = counter + 1
    } else {
        print("succeeded")
        break;
    }
}
##browser()
if(class(x)=="try-error")
{
    print(inFile)
    ## print(x)
    ## print(dir(dirname(inFile)))
    stop("failed to load input")
}


for(funcFileElem in funcFile)
{
    source(funcFileElem)
}

for (aname in c("outFile", "funcArgs", "func", "funcFile"))
{ 
    print("*********")
    print(aname)
 ##   print(head(eval(parse(text=aname))))
}

load(funcFileForObj)
load(otherGlobalsFile)


clusterOut = do.call(func, c(funcArgs, otherGlobals))

success = F
slepe = 5
trials = 3
count = 1
while(!success & count<trials)
{
    x = try(save(clusterOut, file=outFile))
    if(class(x)!="try-error")
    {
        success = T
        break;
    }
    count = count + 1
    slepetime = as.integer(slepe + ceil(abs(rnorm(1, mean = 0, sd = slepe/4))))
    print(paste0("sleeping: ",slepetime))
    Sys.sleep(slepetime)
}

if(!success)
{
    save(clusterOut, file=outFile)
}
