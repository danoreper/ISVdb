library(stringr)
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

writeParams <- function(propObj, filename)
{
    astring = yaml::as.yaml(as.list(propObj))
    write(astring, file = filename)
}


## loads params from command line
loadParams <-  function(args = commandArgs(T))
{

    yamldefault  = "../config/defaultParams.yaml" 
    prop          = yaml::yaml.load_file(yamldefault)

    argz = parseCommandArgs(args)

##    browser()
    
    print(paste0("override arg is ",argz$propArgs))
    if(length(args)>0)
    {
        for(arg in argz$propArgs)
        {
            overrideYaml = arg
            print(arg)
            propOver = yaml::yaml.load_file(overrideYaml)
            prop = .loadParamRecurseHelp(prop,propOver)
        }
    }
    prop = list2env(prop)
    return(prop)
}

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
