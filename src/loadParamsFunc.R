library(stringr)

##Parse command line args; anything without a -I is considred to be a property override
parseCommandArgs <- function(args = commandArgs(T))
{
   ##browser()
    inpArgs  = grepl(args, pattern="^[-I]")
    propArgs = !inpArgs

    inpArgs = args[inpArgs]
    inpArgs = gsub(inpArgs, pattern = "-I", replacement = "")
    inpArgs = unlist(inpArgs)


    propArgs = args[propArgs]
    ## print("parsing!")
    ## print(inpArgs)
    ## print(propArgs)
    return(list(inpArgs = inpArgs, propArgs = propArgs))
}

##write yaml properties to file
writeParams <- function(propObj, filename)
{
    astring = yaml::as.yaml(as.list(propObj))
    write(astring, file = filename)
}


## function for loading params from command line args; the files
loadParams <-  function(args = commandArgs(T))
{

    yamldefault  = "../config/defaultParams.yaml"
    if(!file.exists(yamldefault))
    {
        print("running without properties")
        prop = list()
        return(prop)
    }
    
    prop          = yaml::yaml.load_file(yamldefault)

    argz = parseCommandArgs(args)

##    browser()
    
##    print(paste0("override arg is ",argz$propArgs))
    if(length(args)>0)
    {
        for(arg in argz$propArgs)
        {
            overrideYaml = arg
##            print(arg)
            propOver = yaml::yaml.load_file(overrideYaml)
            prop = .loadParamRecurseHelp(prop,propOver)
        }
    }
    prop = list2env(prop)
    return(prop)
}

##recursive helper for overriding properties; i.e. overriding sublists
.loadParamRecurseHelp <- function(origlist, overlist)
{
    for(aname in ls(overlist))
    {
        if(aname %in% ls(origlist))
        {
            if(!(is.list(origlist[[aname]])))
            {
                origlist[[aname]] = overlist[[aname]]
            } else {
                origlist[[aname]] = .loadParamRecurseHelp(origlist[[aname]], overlist[[aname]])
            }
        }
        else
        {
            origlist[[aname]] = overlist[[aname]]
            warning(paste0("prop in override file doesn't exist in default: ", substitute(overlist),":", aname))
        }
    }
    return(origlist)
}
