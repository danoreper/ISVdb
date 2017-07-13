library(parallel)
source("./utils.R")

##TODO: add the option of deleting temp folders which are successfully completed.
##TODO: consider writing special code for saved data frames.
parallel = new.env(hash=T)

unlink(".RData")

parallel$get.cluster.accum <- function(system.type,       
                                       func,
                                       sharedVariables        = list(),
                                       
                                       filesToSource          = c(),
                                       batchSize              = 1,
                                       timeLimit.hours        = 1,
                                       cpuMemLimit.GB         = 4,
                                       coresPerJob            = 1,
                                       maxSimulJobs           = 1000,
                                       systemOpts             = c(),
                                       outdir                 = fp("./job.tmp"),
                                       retryFailing           = F,
                                       saveProp               = F)
{
    ##Consider getting rid of these as arguements to change.... probably a good idea to package them up into a class like thing, including the script.

    getInFile     = parallel$.getDefaultInputFile
    getOutFile    = parallel$.getDefaultOutputFile
    getFailFile   = parallel$.getDefaultFailFile
    getROutFile   = parallel$.getDefaultROutFile
    getPropFile   = parallel$.getDefaultPropFile
    clusterScript = "./parallel/jobscript.R"
    accum = new.env(hash=T)

    accum$.jobSystem = eval( parse(text = system.type))
    accum$.jobSubmitCommand = accum$.jobSystem$getSubmitCommand(time.hours = timeLimit.hours,
                                                                numProcessPerNode = coresPerJob,
                                                                memoryLimit.GB = cpuMemLimit.GB)
    if(!is.null(systemOpts)) {accum$.jobSubmitCommand = paste(accum$.jobSubmitCommand, paste(systemOpts, collapse = " "))}   
    
    accum$.sleepCheck=15
    accum$.filesToSource = filesToSource
    accum$.batchSize = batchSize
    accum$.jobsLimit = maxSimulJobs

  

    accum$.currentBatchArgs = list()
    accum$.clusterCommands = c()
    funcname = as.character(substitute(func))
    funcname = funcname[length(funcname)]
    x = fp(outdir, paste0(tempdir(), "_",funcname ,"_", parallel$.get.start()))
    accum$.outdir = x

    dir.create(accum$.outdir, showWarnings = F, recursive=T)

    accum$.batchLengths = c()
    accum$.funcFile = fp(accum$.outdir, "func.RData")
    otherGlobals = sharedVariables
    save(list=c("func","otherGlobals"), file = accum$.funcFile)
    
    accum$.ready =  T
    
    accum$addCall <- function(funcArgs, propObj=NULL)
    {
        if(!accum$.ready)
        {
            stop("create a new accumulator, this one has already been run.")
        }
        if(is.null(propObj))
        {
            if(saveProp)
            {
                accum$.propObj = prop ##use the default property object defined globally.
            }
        }
        
        ##IF the current Batch is full, save it, start a new batch
        accum$.currentBatchArgs = util$appendToList(accum$.currentBatchArgs, funcArgs)
        if(length(accum$.currentBatchArgs)==(accum$.batchSize))
        {
            accum$.save.current.batch()
            
        } 
    }
    
    accum$.save.current.batch <- function()
    {
        
        if(length(accum$.currentBatchArgs)==0)
        {
            return()
        }
        i = length(accum$.clusterCommands)+1
        commandRoot = paste0("R CMD BATCH --no-save --no-restore '--args ")
        
        propFile = getPropFile(accum$.outdir, i)
        if(saveProp)
        {
            writeParams(accum$.propObj, filename = propFile)
        }
        
        inFile   = getInFile  (accum$.outdir, i)
        outFile  = getOutFile (accum$.outdir, i)
        failFile = getFailFile(accum$.outdir,i)
        routFile = getROutFile(accum$.outdir,i)
        
        funcArgs = accum$.currentBatchArgs
        funcFile = accum$.funcFile
        filesToSource = accum$.filesToSource
        force(funcArgs)

        ##TODO: change the script to use -P in front of the prop file? or perhaps encode the propFile into the infile, rather than the prop? Or remove the prop?
        save(file = inFile, list = c("outFile","failFile", "funcArgs",  "funcFile", "prop", "filesToSource"))
        command = paste0(commandRoot, " ", propFile, " -I",inFile, "' ",clusterScript, " ", routFile)
##        print(command)
        accum$.clusterCommands <<- c(accum$.clusterCommands, command)
        accum$.batchLengths = c(accum$.batchLengths, length(accum$.currentBatchArgs))
        accum$.currentBatchArgs <<-list()
    }
    
    accum$runAll <- function()
    {
        print(paste0("running all jobs; temp files for debugging stored in ", accum$.outdir))
        accum$.save.current.batch()
        
        localf = function(clusterCommands)
        {
  
            clusterCommands =
                parallel$.submitCommands(jobSystem        = accum$.jobSystem,
                                         jobSubmitCommand = accum$.jobSubmitCommand,
                                         clusterCommands  = clusterCommands,
                                         sleepCheck       = accum$.sleepCheck,
                                         cluster.outdir   = fp(accum$.outdir, "jobsubmit.outfiles"),
                                         outdir           = fp(accum$.outdir),
                                         batchSize        = accum$.batchSize,
                                         maxNumJobs       = accum$.jobsLimit)$failingcommands
        }

        prevFailingCount = Inf
        failing = localf(clusterCommands = accum$.clusterCommands)
        if(retryFailing)
        {
            while(length(failing)>0 && length(failing)<prevFailingCount)
            {
                print(paste0("failures:", length(failing)))
                prevFailingCount = length(failing)
                failing = localf(failing)
            }
        }
        
        outfiles = getOutFile(accum$.outdir, 1:length(accum$.clusterCommands)) 
        accum$.ready <<- F

        return(outfiles)
    }

    accum$getOutputIterator <- function(outputs)
    {
        i  = 1
        j  = 1
        batch.i = parallel$.getOutput(outputs[[i]])
        batchlen.i = accum$.batchLengths[i]

        iter = new.env(hash=T)

        iter$nextItem <- function()
        {
            out     = batch.i[[j]]  
            j <<- j + 1
            if(j>batchlen.i)
            {
                i <<- i + 1
                j <<- 1
                
                if(i>length(outputs))
                {
                    batch.i<<-NULL
                    j <<- 0
                    return(out)
                }

                batch.i <<- parallel$.getOutput(outputs[[i]])
                batchlen.i <<- accum$.batchLengths[i]
            }
            return(out)
        }
        
        iter$hasNext <- function()
        {
            return(!is.null(batch.i)&& !length(batch.i)==0)
        }
        
        return(iter)
    }
    
    
##TODO implement the iterator functionality
    accum$getAllOutputs <- function(outs, removeFailing = F)
    {
        counter = 1
        outputs = list()
        for(i in 1:length(outs))
        {
            outputbatch  = try(parallel$.getOutput(outs[[i]]))
            
            if(is.null(outputbatch)||class(outputbatch)=="try-error")
            {
                print(paste0("failing batch: ",outs[[i]]))
                if(!removeFailing)
                {
                    for(j in 1:accum$.batchLengths[i])
                    {
                        outputs[[counter]] = NA
                        counter = counter + 1
                    }
                }
            } else {
                for(j in 1:length(outputbatch))
                {
                    output = outputbatch[[j]]

                    validout = T
                    if(length(output)==1)
                    {
                        if(is.na(output) || class(output)=="try-error")
                        {
                            validout = F
                        } 
                    }
                    if(!removeFailing || validout)
                    {
                        outputs[[counter]] = output
                        counter = counter + 1
                    }
                }
            }
        }
        return(outputs)
    }
    
    return(accum)
}

    
parallel$.getOutput <- function(outfile)
{
    load(outfile)
    ##This ('clusterOut') is the name of the variable holding the result as defined in generated bsubScript
    return(clusterOut)
}


##TODO: move mclapply stuff to another file
parallel$get.mc.accum <- function(func, mc.cores, sharedVariables = list(), batchSize=100*mc.cores)
{
    accum = new.env(hash=T)
    accum$ready = F
    
    accum$func=func
    accum$sharedVariables = sharedVariables
    accum$funcArgs = list()
    accum$ready = T

    accum$.funcWrapper <- function(ind)
    {
        argz = accum$funcArgs[[ind]]
        argz = c(argz, accum$sharedVariables)
        out = try(do.call(accum$func, argz))

        return(out)
    }
        
    accum$addCall <- function(funcArgs)
    {
        if(!accum$ready)
        {
            stop("init should have beein called after results were collected, or at the beginning before any calls")
        }
        ##        i = length(accum$funcArgs)+1
        accum$funcArgs <<- util$appendToList(accum$funcArgs, funcArgs)
    }

    accum$runAll <- function()
    {
        inds = 1:length(accum$funcArgs)
        if(length(inds)==1 & class(inds)=="logical")
        {
            browser()
        }
        out = parallel$lapply.wrapper(inds, FUN = accum$.funcWrapper, batchSize = batchSize, mc.cores = mc.cores)
        print("ranall")
        accum$ready = F
        return(out)
    }

    accum$getAllOutputs <- function(outs, removeFailing = F)
    {
        if(!removeFailing)
        {
            return(outs)
        } else {

            cleaned = list()
            counter = 1
            for(elem in outs)
            {
                if((length(elem)==1 && is.na(elem))||class(elem)=="try-error")
                {
                    next
                } else {
                    cleaned[[counter]] = elem
                    counter = counter + 1
                }
            }
            return(cleaned)
        }
    }

    accum$getOutputIterator <- function(outputs)
    {
        i  = 1
        iter = new.env(hash=T)

        iter$nextItem <- function()
        {
            if(i>length(outputs))
            {
                out = NULL
            } else {
                out = outputs[[i]]
                i <<- i+1
            }
            return(out)
        }

        iter$hasNext <- function()
        {
            return(i<=length(outputs))
        }
        
        return(iter)
    }

    return(accum)
}

 
##TODO move lapply stuff to a separate file
##FUN must take as its first argument a vector of indexes
##and grab the relevant portion of whatever the additional arguments are
parallel$lapply.wrapper <- function(inds, FUN, batchSize = 10*mc.cores, mc.cores, ...)
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



##Need to revamp the logging to be more useful

##submit a list of commands prefixed by the bsubCommand, and if sleepcheck is not null,
##blocks untill all commands complete, checking at an interval of sleepCheck seconds to
##see if this is so. 
##TODO: incorporate max jobs
parallel$.submitCommands <- function(jobSystem, jobSubmitCommand, clusterCommands, sleepCheck=30,  cluster.outdir, outdir,  batchSize = 1, maxNumJobs)
{
    pracma::tic()
    start.time  = parallel$.get.start()
    ## a map from jobname, to the command or commands that are called by that job
    submitted.jobs = list()
    for(i in 1:length(clusterCommands))
    {
        while(jobSystem$countActiveJobs(names(submitted.jobs))>=maxNumJobs)
        {
            numActive = jobSystem$countActiveJobs(names(submitted.jobs))
            print(paste0("sleeping until fewer than ", maxNumJobs, " are active; currently active=", numActive))
            Sys.sleep(sleepCheck)
        }
        
        jobname        = paste0(start.time, ".", i)
        submitted.jobs[[jobname]] = clusterCommands[i]

        jobSystem$run.single.job(jobSubmitCommand, jobname, clusterCommands[i], cluster.outdir)
    }
    print(paste0("submitted all jobs started at ", start.time))

    if(!is.null(sleepCheck))
    {
        while(jobSystem$countActiveJobs(names(submitted.jobs))>0)
        {
            numActive = jobSystem$countActiveJobs(names(submitted.jobs))
            print(paste0("sleeping untill the remaining ", numActive ," jobs are finished"))
            Sys.sleep(sleepCheck) 
        }
    }
    print("done with jobs")
    pracma::toc()

    
    failingcommands = c()
    failingjobs = c()
    for (i in 1:length(submitted.jobs))
    {
        jobname = names(submitted.jobs)[[i]]
        outfile = parallel$.getOutLogFile(cluster.outdir,jobname)

        failedBatch =
            (!file.exists(outfile) && jobSystem$thename == "killdevil") ||
            (!file.exists(parallel$.getDefaultOutputFile(cluster.outdir, jobname)) && jobSystem$thename == "longleaf")
        
        
        if(failedBatch)
        {
            print("********************")
            print("Failed entire batch, no outfile:")
            print(jobname)
            print(submitted.jobs[[jobname]])
            print("********************")

            failingcommands = c(failingcommands, submitted.jobs[[jobname]])
            failingjobs     = c(failingjobs, jobname)
            next
        }

        if(jobSystem$thename == "killdevil")
        {
            success = jobSystem$checkCompletion(outfile)
        } else {
            success  = jobSystem$checkCompletion(jobname)
        }
    
        if(!success)
        {
            failingcommands = c(failingcommands, submitted.jobs[[jobname]])
            failingjobs     = c(failingjobs, jobname)
            print("********************")
            print("Failed entire batch:")
            print(jobname)
            print(submitted.jobs[[jobname]])
            print("********************")
            next
        }

        failfile = parallel$.getDefaultFailFile(outdir, i)
        load(failfile)
        if(length(failed)>1)
        {
            print("********************")
            print("Failed runs within batch:")
            print(jobname)
            print(submitted.jobs[[jobname]])
            print(failed)
            print("********************")
            failingcommands = c(failingcommands, submitted.jobs[[jobname]])
            failingjobs     = c(failingjobs, jobname)
        }
              
    }
    
    return(list(failingcommands = failingcommands, failingjobs = failingjobs))
}




##TODO: make these private?
parallel$.get.start <- function()
{
    start.time = gsub(pattern = " ", replacement  ="_",format(Sys.time()))
    return(start.time)
}

parallel$.getOutLogFile <- function(outputlocaldir,jobname)
{
    return(fp(outputlocaldir, paste0(jobname, ".job.out")))
}

parallel$.getErrorLogFile <- function(outputlocaldir,jobname)
{
    return(fp(outputlocaldir, paste0(jobname, ".job.err")))
}


##input file contain the function arguments
parallel$.getDefaultInputFile <- function(outdir, i)
{
    outdirin = fp(outdir, "in")
    dir.create(outdirin, showWarnings = F, recursive = T)
    return(fp(outdirin, paste0("in_", i)))
}

##prop files contain the property overrides specific to this run 
parallel$.getDefaultPropFile <- function(outdir, i)
{
    outdirin = fp(outdir, "prop")
    dir.create(outdirin, showWarnings = F, recursive = T)
    return(fp(outdirin, paste0("prop_", i)))
}

##output files contain the output of calling the function
parallel$.getDefaultOutputFile <- function(outdir, i)
{
    outdirout = fp(outdir, "out")
    dir.create(outdirout, showWarnings = F, recursive = T)
    return(fp(outdirout, paste0("out_", i)))
}

##output files contain the output of calling the function
parallel$.getDefaultFailFile <- function(outdir, i)
{
    outdirout = fp(outdir, "fail")
    dir.create(outdirout, showWarnings = F, recursive = T)
    return(fp(outdirout, paste0("fail_", i)))
}


##Rout files contain the .ROut files specific to a function call.
parallel$.getDefaultROutFile <- function(outdir, i)
{
    outdirout = fp(outdir, "ROut")
    dir.create(outdirout, showWarnings = F, recursive = T)
    return(fp(outdirout, paste0("out_", i, ".ROut")))
}

parallel$.run.single.job <- function(submitCommand, jobname, quantcommandset, outputlocaldir, outflag, errorflag, quoteCommand)
{
    submitCommand = paste0(submitCommand, " -J ", jobname)
    if(!is.null(outputlocaldir))
    {
        dir.create(outputlocaldir, showWarnings=F, recursive =T)
        submitCommand = paste(submitCommand, outflag, " ", parallel$.getOutLogFile(outputlocaldir, jobname))
        submitCommand = paste(submitCommand, errorflag, " ", parallel$.getErrorLogFile(outputlocaldir, jobname))
    }
    if(quoteCommand)
    {
        qbreak = " \" "
    } else { qbreak = ""}
    fullcommand = paste(submitCommand, qbreak,  paste(quantcommandset, collapse="; "), qbreak)
##    cat(fullcommand)
    x = invisible(system(fullcommand,intern=T, ignore.stdout = T))
    return(x)
}




killdevil = new.env(hash = T)
killdevil$thename = "killdevil"

killdevil$getSubmitCommand <- function(time.hours, numProcessPerNode, memoryLimit.GB)
{
    if(time.hours<1)
    {
        queue = "hour"
    } else if (time.hours<24) {
        queue = "day"
    } else {
        queue = "week"
    }
        
    command =  paste0("bsub -R 'span[hosts=1]' -n ", numProcessPerNode, " -M ", ceiling(memoryLimit.GB) ," -q ", queue)
}

    
##Needs to change for longleaf
killdevil$run.single.job <- function(submitCommand, jobname, quantcommandset, outputlocaldir=NULL)
{
    command = parallel$.run.single.job(submitCommand, jobname, quantcommandset, outputlocaldir, outflag="-oo", errorflag="-eo", quoteCommand = T)
    return(command)
}


killdevil$countActiveJobs <- function(submitted.jobs)
{
    if(is.null(submitted.jobs))
    {
        return(0)
    }
    a = try(system("bjobs -w", intern=T))
    
    tokens = strsplit(a[1], "\\s+")[[1]]
    colind = which(tokens == "JOB_NAME") 
    jobids = strsplit(a[2:length(a)], "\\s+")
    jobids = unlist(lapply(jobids, "[", colind))
    numActive = length(intersect(jobids, submitted.jobs))
    return(numActive)
}

killdevil$checkCompletion <- function(outfile)
{
    grepcom = paste("grep -l ", "'Successfully completed.'", outfile)
    out = try(system(grepcom, intern=T))
    if(class(out)=="try-error"||length(out)==0)
    {
        return(F)
    } else {
        return(T)
    }
}



longleaf = new.env(hash = T)
longleaf$thename = "longleaf"

longleaf$getSubmitCommand <- function(time.hours, numProcessPerNode, memoryLimit.GB)
{

    time.minutes = ceiling(time.hours*60)
    mem.MB = ceiling(memoryLimit.GB*1024)
    
    command = paste0("sbatch -t ", time.minutes, " -n ", numProcessPerNode, " --mem-per-cpu=",mem.MB)
}

##package the three methods below?
##Needs to change for longleaf
longleaf$run.single.job <- function(submitCommand, jobname, quantcommandset, outputlocaldir=NULL)
{
    command = parallel$.run.single.job(submitCommand, jobname, quantcommandset, outputlocaldir, outflag="-o", errorflag="-e", quoteCommand = F)
}

longleaf$checkCompletion <- function(jobname)
{
    theuser = system("echo $USER", intern = T)
    statuscommand = paste0("sacct -n --format State -u ", theuser, " --name=", jobname)
  
    out = try(system(statuscommand, intern=T))
    if(class(out)=="try-error"||length(out)==0)
    {
        return(F)
    } else {
        if(all(out== " COMPLETED "))
        {
            return(T)
        } else {
            print("failed job in check completion!!!")
            print(out)
            return(F)
        }
    }
}

longleaf$countActiveJobs <- function(submitted.jobs)
{
    theuser = system("echo $USER", intern = T)
    command = paste0("squeue -u ", theuser, " -h ", " --name=", paste(submitted.jobs, collapse=","), "|wc -l")
    numJobs = as.numeric(system(command, intern = T))
    return(numJobs)
}


