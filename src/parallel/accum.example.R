source("./loadParams.R")
source("./parallel/accumulator.R")
##source("./utils.R")



f = function(x,y)
{
    Sys.sleep(.1)
    if(x<10)
    {
        return(asdfafds)
    }
    return(x+y)
}

tmpdir = "./job.tmp/"

##system.type = "longleaf"
system.type = "mc"
##system.type = "killdevil"

if(system.type == "mc")
 {
     accum = parallel$get.mc.accum(func = f, mc.cores = 1, sharedVariables = list(y=100))
 } else {

     accum = parallel$get.cluster.accum(system.type       = system.type,
                                        func = f,
                                        sharedVariables   = list(y= 100),
                                        
                                        filesToSource     = c(),
                                        batchSize         = 7,
                                        timeLimit.hours   = .15,
                                        cpuMemLimit.GB    = .5,
                                        coresPerJob       = 1,
                                        ##                                   maxSimulJobs      = 4,
                                        systemOpts        = c(),
                                        outdir            = tmpdir,
                                        retryFailing      = T,
                                        saveProp          = T)
 }


for(i in 1:30)
{
    accum$addCall(list(x = i))
}

outputs = accum$runAll()

browser()


print("getAllOutputs")
print(unlist(accum$getAllOutputs(outputs)))

print("getAllOutputs removing failing")
print(unlist(accum$getAllOutputs(outputs, removeFailing = T)))


print("iterate over outputs")
iter = accum$getOutputIterator(outputs)
k = 1
## debug(iter$hasNext)
## debug(iter$nextItem)

while(iter$hasNext())
{
    item = iter$nextItem()
    print(item)
}


## SLURM
## library(rslurm)

## f = function(x)
## {
    
##     z = function(x)
##     {
##         return(x*x)
##     }
    
##     return(z(x))
## }


## k = slurm_call(f, params = list(x=5))
## out = get_slurm_out(k)
## print(out)


