##
## Functions for parallelizing both on the same computer, and for submitting jobs via bsub to killdevil cluster
##
#################################################################
library(parallel)
source("./utils.R")

bsub = new.env(hash=T)


tmpDir = fp(prop$tmpdir)
bsub$get.bsub.accumulator <- function(funcFile, bsubCommand,
                                      batchSize     = 1,
                                      baseOut       = tmpDir,
                                      getInFile     = bsub$getDefaultInputFile,
                                      getOutFile    = bsub$getDefaultOutputFile,
                                      getROutFile   = bsub$getDefaultROutFile,
                                      getPropFile   = bsub$getDefaultPropFile,
                                      clusterScript = "./bsubScript.R")
{
    sleepCheck=30
    accumulator = new.env(hash=T)
    accumulator$ready = F
    
    accumulator$init <- function(func, otherGlobals=list())
    { 

        accumulator$quantCommands <<- c()
##        try(unlink(outdir, recursive = T), showWarnings = F)
        outdir                    <<- fp(baseOut, paste0(tempdir(), bsub$get.start()))
        dir.create(outdir, showWarnings = F, recursive=T)
        accumulator$ready         <<- T
        funcFileForObj <<- fp(outdir, "func.RData")
        unlink(funcFileForObj)
        save(func, file = funcFileForObj)
        
        otherGlobalsFile <<- fp(outdir, "otherGlobals.RData")
        unlink(otherGlobalsFile)
        save(otherGlobals, file = otherGlobalsFile)
    }
    
    accumulator$addCall <- function(funcArgs, propObj=NULL)
    {
        if(!accumulator$ready)
        {
            stop("init should have beein called after results were collected, or at the beginning before any calls")
        }
        ## srcfile = attr(attr(func, "srcref"), "srcfile")$filename
        ## print(srcfile)
        if(is.null(propObj))
        {
            propObj = prop ##use the default property object defined globally.

        }
                
        commandRoot = paste0("R CMD BATCH --no-save --no-restore '--args ")
        i = length(accumulator$quantCommands)+1
        propFile = getPropFile(outdir, i)
        writeParams(propObj, filename = propFile)
        inFile   = getInFile  (outdir, i)
        outFile  = getOutFile (outdir, i)
        routFile = getROutFile(outdir,i)
##        funcArgs = eval(parse(text="funcArgs"))
        save(file = inFile, list = c("outFile", "funcArgs",  "funcFile", "prop", "funcFileForObj", "otherGlobalsFile"))

        command = paste0(commandRoot, " ", propFile, " -I",inFile, "' ",clusterScript, " ", routFile)
##        print(command)
        accumulator$quantCommands <<- c(accumulator$quantCommands, command)
    }

    accumulator$runAll <- function()
    {
        logdir = fp(outdir, "bsub.outfiles")

        prevFailingCount = Inf
        failing = bsub$submitCommands(bsubCommand, accumulator$quantCommands, sleepCheck, logdir, batchSize=batchSize)$failingcommands

        while(length(failing)>0 & length(failing)<prevFailingCount)
        {
            print(failing)
            print(paste0("failures:", length(failing)))
            prevFailingCount = length(failing)
            failing = bsub$submitCommands(bsubCommand, failing, sleepCheck, logdir, batchSize=batchSize)$failingcommands
        }

        
        outfiles = getOutFile(outdir, 1:length(accumulator$quantCommands)) 
        accumulator$ready <<- F
        
        return(outfiles)
    }

    accumulator$getOut <- function(outobj)
    {
        bsub$getOutput(outobj)
    }
    
    accumulator$outputs.files = T

    return(accumulator)
}


bsub$getAllOutputs <- function(outs, accum)
{
    outputs = list()
    for(i in 1:length(outs))
    {
        output  = try(accum$getOut(outs[[i]]))

        if(is.null(output)||class(output)=="try-error")
        {
            print(paste0("failing outfile: ",outs[[i]]))
            outputs[[i]] = NA
        } else {
            outputs[[i]] = output
        }
            
    }
    return(outputs)
}

bsub$getOutput <- function(outfile)
{
    load(outfile)
    ##This ('clusterOut') is the name of the variable holding the result as defined in generated bsubScript
    return(clusterOut)
}

bsub$get.stub.accum <- function()
{
    return(bsub$get.mc.accumulator(1))
}

bsub$get.mc.accumulator <- function(mc.cores, batchSize=100*mc.cores)
{
    accumulator = new.env(hash=T)
    accumulator$ready = F
    
    accumulator$.funcWrapper <- function(ind)
    {
##        print("calling func wrapper")
        ## out = list()
        ## for(i in inds)
##        {
        argz = accumulator$funcArgs[[ind]]
        argz = c(argz, accumulator$otherGlobals)
        out = do.call(accumulator$func, argz)
        ##      }
        return(out)
    }

    accumulator$init <- function(func, otherGlobals = list())
    {
        accumulator$func<<-func
        accumulator$otherGlobals <<- otherGlobals
        accumulator$funcArgs <<- list()
        accumulator$ready = T
    }
    
    accumulator$addCall <- function(funcArgs)
    {
        if(!accumulator$ready)
        {
            stop("init should have beein called after results were collected, or at the beginning before any calls")
        }
        ##        i = length(accumulator$funcArgs)+1
        accumulator$funcArgs <<- util$appendToList(accumulator$funcArgs, funcArgs)
    }

    accumulator$runAll <- function()
    {
        inds = 1:length(accumulator$funcArgs)
        if(length(inds)==1 & class(inds)=="logical")
        {
            browser()
        }
        out = bsub$lapply.wrapper(inds, FUN = accumulator$.funcWrapper, batchSize = batchSize, mc.cores = mc.cores)
##        print("ranall")
        accumulator$ready = F
        return(out)
    }

    accumulator$getOut <- function(out)
    {
        return(out)
    }
    
    accumulator$outputs.files = F

    return(accumulator)
    
}

 
##TODO move to parallel env
##FUN must take as its first argument a vector of indexes
##and grab the relevant portion of whatever the additional arguments are
bsub$lapply.wrapper <- function(inds, FUN, batchSize = 10*mc.cores, mc.cores, ...)
{
    
    if(mc.cores==1)
    {
        results = lapply(X=inds, FUN=FUN, ...)
        return(results)
    } 
    else 
    {
        ##Split into batches to allow mclapply to work properly, it doesn't garbage collect well.
        ##function call returns a list of lists of indexes, splitting up inds into indexes of length batch size
        indexGroups = util$getIndexGroupsForInds(inds, batchSize) 

        results = list()
        for(i in 1:length(indexGroups))
        {
            indsGroup = indexGroups[[i]]
            print(paste0("working on index: ", indsGroup[1]))
            results[[i]] = mclapply(X=indsGroup, FUN=FUN, mc.cores = mc.cores, ...) 
        }

        ## at this point we have a list of lists (the outermost list corresponds to separate batches) and we'd like the returned value to avoid reflecting the innards of this method-- merge all the batches together

        results = do.call(c, results)
        return(results)
    }
    ##eliparg = list(...)
}



##submit a list of commands prefixed by the bsubCommand, and if sleepcheck is not null,
##blocks untill all commands complete, checking at an interval of sleepCheck seconds to
##see if this is so. Sleeps rest of the time, so probably ok to run on the cluster
bsub$submitCommands <- function(bsubCommand, quantCommands, sleepCheck=NULL, bsuboutdir = NULL, batchSize = 1)
{
    len = length(quantCommands)
    indexGroups = util$getIndexGroupsForLen(len, batchSize)

    start.time  = bsub$get.start()
    ## a map from jobname, to the command or commands that are called by that job
    submitted.jobs = list()
    for(i in 1:length(indexGroups))
    {
        indsGroup = indexGroups[[i]]
        quantCommandSet = quantCommands[unlist(indsGroup)]
        jobname        = paste0(start.time, ".", i)
        submitted.jobs[[jobname]] = quantCommandSet
        bsub$run.single.bsub(bsubCommand, jobname, quantCommandSet,bsuboutdir)
    }

    
    if(!is.null(sleepCheck))
    {
        bsub$block.on.bsub(names(submitted.jobs), sleepCheck)
    }

    failingcommands = c()
    failingjobs     = c()
    if(!is.null(bsuboutdir))
    {
        for (jobname in names(submitted.jobs))
        {
            outfile = bsub$getOutLogFile(bsuboutdir,jobname)
            if(!file.exists(outfile))
            {
                print("job failed!")
##                print(submitted.jobs[[jobname]])
                
                failingcommands = c(failingcommands, submitted.jobs[[jobname]])
                failingjobs     = c(failingjobs, jobname)
                next
            }

            grepcom = paste("grep -l ", "'Successfully completed.'", outfile)
            out = try(system(grepcom, intern=T))
            if(class(out)=="try-error")
            {
                failingjobs = c(failingjobs, jobname)
                print(out)
                browser()
                
            } else {
                if(length(out) == 0)
                {
                    failingcommands = c(failingcommands, submitted.jobs[[jobname]])
                    failingjobs     = c(failingjobs, jobname)
                }
            }
        }
    }
    return(list(failingcommands = failingcommands, failingjobs = failingjobs))
}


bsub$get.start <- function()
{
    start.time = gsub(pattern = " ", replacement  ="_",format(Sys.time()))
    return(start.time)
}

bsub$getOutLogFile <- function(outputlocaldir,jobname)
{
    return(fp(outputlocaldir, paste0(jobname, ".bsub.out")))
}

bsub$getErrorLogFile <- function(outputlocaldir,jobname)
{
    return(fp(outputlocaldir, paste0(jobname, ".bsub.err")))
}


bsub$run.single.bsub <- function(bsubCommand, jobname, quantcommandset, outputlocaldir=NULL)
{
    bsubCommand = paste0(bsubCommand, " -J ", jobname)
    if(!is.null(outputlocaldir))
    {
        dir.create(outputlocaldir, showWarnings=F, recursive =T)
        bsubCommand = paste(bsubCommand, "-oo ", bsub$getOutLogFile(outputlocaldir, jobname))
        bsubCommand = paste(bsubCommand, "-eo ", bsub$getErrorLogFile(outputlocaldir, jobname))
    }
    fullcommand = paste(bsubCommand,
                        " \"  ",
                        paste(quantcommandset, collapse="; "),
                        " \" ")
    ## cat(fullcommand)
    ## browser()
    system(fullcommand)
}

bsub$block.on.bsub <- function(submitted.jobs, sleepCheck)
{
    while(T)
    {
        Sys.sleep(sleepCheck) 
        a = try(system("bjobs -w", intern=T))
        ##
        if(length(a)==0)
        {
            break;
        }
        tokens = strsplit(a[1], "\\s+")[[1]]
        colind = which(tokens == "JOB_NAME") 
        jobids = strsplit(a[2:length(a)], "\\s+")
        jobids = unlist(lapply(jobids, "[", colind))
        if(length(intersect(jobids, submitted.jobs))==0)
        {
            break;
        } 
    }
}

##We want to make sure the input files are in a separate directory from where the output files are being generated; there may be some problems with multiple nodes reading from and writing to the same directory. The input files are generated serially before any output files are generated

##input file contain the function arguments
bsub$getDefaultInputFile <- function(outdir, i)
{
    outdirin = fp(outdir, "in")
    dir.create(outdirin, showWarnings = F, recursive = T)
    return(fp(outdirin, paste0("in_", i)))
}

##prop files contain the property overrides specific to this run 
bsub$getDefaultPropFile <- function(outdir, i)
{
    outdirin = fp(outdir, "prop")
    dir.create(outdirin, showWarnings = F, recursive = T)
    return(fp(outdirin, paste0("prop_", i)))
}

##output files contain the output of calling the function
bsub$getDefaultOutputFile <- function(outdir, i)
{
    outdirout = fp(outdir, "out")
    dir.create(outdirout, showWarnings = F, recursive = T)
    return(fp(outdirout, paste0("out_", i)))
}

##Rout files contain the .ROut files specific to a function call.
bsub$getDefaultROutFile <- function(outdir, i)
{
    outdirout = fp(outdir, "ROut")
    dir.create(outdirout, showWarnings = F, recursive = T)
    return(fp(outdirout, paste0("out_", i, ".ROut")))
}


bsub$get.default.killdevil.bsub <- function(numProcessPerNode, memoryLimit.GB, queue)
{
    command =  paste0("bsub -R 'span[hosts=1]' -n ", numProcessPerNode, " -M ", memoryLimit.GB ," -q ", queue)
}
