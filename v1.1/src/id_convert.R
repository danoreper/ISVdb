source("./loadParams.R")
id_convert = new.env(hash=T)

id_convert$old.cc.to.new <- function(old.cc)
{
    mapping = fread(dat(prop$genome$ccNameMap))
    newids = rep(NA, length(old.cc))
    for(i in 1:length(old.cc))
    {
        old.cc.elem = old.cc[i]
        inds = grepl(old.cc.elem, mapping$Alias)
        if(sum(inds)!=1)
        {
            stop(paste0("less or more than 1 id corresponding to ", old.cc.elem, ":", paste(mapping$alias[inds], collapse=",")))
        }
        newids[i] = mapping$StrainName[inds]
        newids[i] = unlist(lapply(strsplit(newids[i], "/") ,"[", 1))
    }

    ## grepl(mapping$Strain.name, old.cc)
    ## setkey(mapping,"Alias")
    ## new.cc = mapping[old.cc]$StrainName
    ## return(new.cc)
    ## 
    return(newids)
}

id_convert$new.cc.to.old <- function(new.cc)
{
    mapping = fread(dat(prop$genome$ccNameMap))
    setkey(mapping,"StrainName")
    new.cc = mapping[new.cc]$Alias
    return(new.cc)
}

