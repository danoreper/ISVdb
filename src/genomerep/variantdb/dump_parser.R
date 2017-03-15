source("./bsub.R")
source("./loadParams.R")

dump_parser = new.env(hash=T)

##TODO: property?
dump_parser$strainsdir = dat(fp("dump_100bp_2017-03-11/genotype"))

##dump_parser$strainsdir = dat(fp("dump_100bp_2017-03-11/diplotype"))

dump_parser$iterate <- function(parseFunc,
                                strains = dump_parser$getStrains(),
                                chrs = dump_parser$getChrs(strain = strains[1]),
                                zipped = T,
                                accum = bsub$get.mc.accumulator(mc.cores = 1),
                                iterator = F)    
{
    outs = list()
    force(parseFunc)


    accum$init(func = dump_parser$callparse, otherGlobals= list(parseFunc = parseFunc, zipped = zipped))
    for(strain in strains)
    {
        for(chr in chrs)
        {
            ##  outs[[paste0(strain,",",chr)]] = parseFunc(strain, chr, df)
            ##            print(paste0(strain,"_",chr))

            force(parseFunc)
            accum$addCall(funcArgs = list(strain = strain, chr = chr))
        }
    }

    outs = accum$runAll()
    
    
    if(iterator)
    {
        stop("unimplemented")
    } else {
        outs = bsub$getAllOutputs(outs, accum)
        for(i in 1:length(outs))
        {
            if(length(outs[[i]])==1 &is.na(outs[[i]]))
            {
                print(paste0("fixing:",i))
                outs[[i]] = NULL
            }
            if(i==length(outs))
            {
                break
            }
        }
        print("alldone!")
    }

    print("loaded up outs")
    return(outs)
}


## dump_parser$iterate.dip <- function(parseFunc,
##                                 strains = dump_parser$getStrains(),
##                                 chrs = dump_parser$getChrs(strain = strains[1]),
##                                 zipped = T,
##                                 accum = bsub$get.mc.accumulator(mc.cores = 1),
##                                 iterator = F)    
## {
##     outs = list()
##     force(parseFunc)


##     accum$init(func = dump_parser$callparse, otherGlobals= list(parseFunc = parseFunc, zipped = zipped))
##     for(strain in strains)
##     {
##         for(chr in chrs)
##         {
##             ##  outs[[paste0(strain,",",chr)]] = parseFunc(strain, chr, df)
##             ##            print(paste0(strain,"_",chr))

##             force(parseFunc)
##             accum$addCall(funcArgs = list(strain = strain, chr = chr))
##         }
##     }

##     outs = accum$runAll()
    
    
##     if(iterator)
##     {
##         stop("unimplemented")
##     } else {
##         outs = bsub$getAllOutputs(outs, accum)
##         for(i in 1:length(outs))
##         {
##             if(length(outs[[i]])==1 &is.na(outs[[i]]))
##             {
##                 print(paste0("fixing:",i))
##                 outs[[i]] = NULL
##             }
##             if(i==length(outs))
##             {
##                 break
##             }
##         }
##         print("alldone!")
##     }

##     print("loaded up outs")
##     return(outs)
## }



dump_parser$getStrains <- function()
{
    founders = c("AJ", "B6", "129", "NOD", "NZO", "CAST", "PWK", "WSB")

    names(founders) = c("A_J", "C57BL6J", "129S1_SvImJ", "NOD_ShiLtJ", "NZO_HlLtJ", "CAST_EiJ", "PWK_PhJ", "WSB_EiJ")    

    cclines = paste0("CC", sprintf(1:1000, fmt = "%003d"))
    names(cclines) = cclines
    all.lines = c(names(founders), cclines)
    
    strains = dir(dump_parser$strainsdir, full.names=F)
    
##    strains = all.lines[strains]
    strains = factor(strains, levels = all.lines)
    strains = sort(strains)
    strains = as.character(strains)
    
##    strains = setdiff(strains, "OR559b38V01")
 ##   strains = c("CC001", "CC002")
    
    return(strains)
}

dump_parser$getChrs <- function(strain = dump_parser$getStrains()[1])
{
    straindir = fp(dump_parser$strainsdir, strain)
    chrfls    = dir(straindir, pattern = ".*\\.txt\\.tar\\.gz")
    chrs      = gsub(chrfls, pattern = "\\.txt\\.tar\\.gz*", replacement = "")
 
    
    chrsAll = c(1:19, "X", "Y", "MT")
    chrs = factor(chrs, levels = chrsAll)
    chrs = as.character(sort(chrs))
    return(chrs)
}


dump_parser$read <- function(strain, chr, zipped=T)
{
    straindir = fp(dump_parser$strainsdir, strain)
    chrfl.long.tgz = fp(straindir, paste0(chr, ".txt.tar.gz"))
    chrfl.long     = gsub(chrfl.long.tgz, pattern = "\\.tar\\.gz", replacement = "")
    if(zipped)
    {
        command = paste0("tar -xOzf ", chrfl.long.tgz)
        df = fread(input=command, stringsAsFactors=T, skip =1)
    } else { 
        df = fread(chrfl.long, stringsAsFactors = T, skip = 1)
    }

    ##TODO remove
    ##  col.names = c("variant_id", "pos", "allele_1", "allele_2", "prob", "is_max", "consequence_1", "consequence_2", "gene_name", "transcript_name")
    ## setnames(df, old = colnames(df), new=col.names)

}



dump_parser$write <- function(strain, chr, df, headerlines)
{
    straindir = fp(dump_parser$strainsdir, strain)
    chrfl.long.tgz = fp(straindir, paste0(chr, ".txt.tar.gz"))
    chrfl.long     = gsub(chrfl.long.tgz, pattern = "\\.tar\\.gz", replacement = "")

    
    cat(headerlines, file=chrfl.long)
    fwrite(df, file = chrfl.long, row.names = F, sep = ",", append =T)

##    system(paste0("rm ", chrfl.long.tgz))
    old.wd = getwd()
    setwd(dirname(chrfl.long))
    
    command = paste("tar -czvf ", basename(chrfl.long.tgz),  basename(chrfl.long))
     
    system(command)
    setwd(old.wd)
    system(paste0("rm ", chrfl.long))
}



dump_parser$callparse <- function(parseFunc, strain, chr, zipped = T)
{

    ##TODO remove
    ## {
    ##     fles = gsub(x=dir(straindir),
    ##             pattern= "_s[0-9]+\\.txt",
    ##             replacement = ".txt")

    ##     file.rename(fp(straindir, dir(straindir)), fp(straindir,fles))

    ##     fles = gsub(x=dir(straindir),
    ##                 pattern= "chr",
    ##                 replacement = "")
    ##     file.rename(fp(straindir, dir(straindir)), fp(straindir,fles))
    ## }

    df = dump_parser$read(strain, chr, zipped = zipped)
    out = parseFunc(strain, chr, df)
    return(out)
}




##                command = paste0("tar -xzvf ", chrfl.long.tgz, " -C ", straindir)
##            print(command)
##                system(command)

    ## if(zipped)
    ## {

    ## }
