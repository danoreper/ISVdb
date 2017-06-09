source("./loadParams.R")
source("./bsub.R")
library(data.table)

db_builder = new.env(hash=T)



db_builder$get.db.lite <- function(dbdir)
{

    dblite = new.env(hash=T)
    dblite$genotype  = db_builder$getInstance(fp(dbdir, "genotype"))
    dblite$diplotype = db_builder$getInstance(fp(dbdir, "diplotype"))
    dblite$genotypeSampling  = db_builder$getInstance(fp(dbdir, "genotypeSampling"))
    dblite$diplotypeSampling = db_builder$getInstance(fp(dbdir, "diplotypeSampling"))

    dblite$temp      = fp(dbdir, "temp") ##db_builder$getInstance(fp(dbdir, "temp") #fp(dbdir, temp)
    dir.create(dblite$temp, showWarnings = F, recursive = T)

    dblite$get.temp.dir <- function(subdir=NULL)
    {
        if(is.null(subdir))
        {
            return(dblite$temp)
        } else {
            adir = fp(dblite$temp, subdir)
            dir.create(adir, recursive=T, showWarnings = F)
            return(adir)
        }
    } 
        
    dblite$iterate <- function(type = "genotype",
                               parseFunc,
                               strain1s= dblite[[type]]$getStrains(), 
                               strain2s= NA,
                               chrs = dblite[[type]]$getChrs(strain1s[1]),
                               select = NULL,
                               zipped = T,
                               accum = bsub$get.mc.accumulator(mc.cores = 1),
                               iterator = F)    
    {
        outs = list()
        force(parseFunc)
        
        accum$init(func = dblite$.callparse, otherGlobals= list(parseFunc = parseFunc, type = type))
        for(i in 1:length(strain1s))
        {
            strain1 = strain1s[i]
            strain2 = strain2s[i]
            for(chr in chrs)
            {
                accum$addCall(funcArgs = list(strain1 = strain1, strain2 = strain2, chr = chr))
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

    dblite$.callparse <- function(parseFunc, strain1, strain2, chr, type)
    {

        if(is.na(strain2))
        {
            inst = dblite[[type]]
            df =  inst$read(strain1, chr)
            out = parseFunc(strain1, chr, df)
        } else {
            inst = dblite[[paste0(type,"Sampling")]]
            df1 =  inst$read(strain1, chr)
            df2 =  inst$read(strain2, chr)

            if(type=="genotype")
            {
                setkey(df1, "variant_id", "transcript_id")
                setkey(df2, "variant_id", "transcript_id")
            } else if (type=="diplotype")
            {
                setkey(df1, "variant_id")
                setkey(df2, "variant_id")
            }
            
            df = df1[df2, prob:=prob*i.prob]
            df[,i.prob:=NULL]
            df[,i.pos:=NULL]
            df[,i.transcript_name:=NULL]
            df[,i.gene_name:=NULL]
            browser()
            out = parseFunc(strain1, strain2, chr, df) 
        }
        return(out)
    }

    return(dblite)
}

db_builder$getInstance <- function(tabledir)
{
    inst.sep = "\t"
    
    dir.create(tabledir, showWarnings = F, recursive = T)
    force(tabledir)
    inst = new.env(hash = T)

    inst$clear <- function()
    {
        unlink(tabledir, recursive = T)
        dir.create(tabledir, showWarnings = F, recursive = T)
    }


    
    
    inst$iterate <- function(parseFunc,
                             strains = inst$getStrains(),
                             chrs = inst$getChrs(strain = strains[1]),
                             select = NULL,
                             zipped = T,
                             accum = bsub$get.mc.accumulator(mc.cores = 1),
                             iterator = F)    
    {
        outs = list()
        force(parseFunc)
        
        accum$init(func = inst$.callparse, otherGlobals= list(parseFunc = parseFunc, zipped = zipped, select = select))
        for(strain in strains)
        {
            for(chr in chrs)
            {
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
    
    inst$getStrains <- function()
    {
        founders = c("AJ", "B6", "129", "NOD", "NZO", "CAST", "PWK", "WSB")
        
        names(founders) = c("A_J", "C57BL6J", "129S1_SvImJ", "NOD_ShiLtJ", "NZO_HlLtJ", "CAST_EiJ", "PWK_PhJ", "WSB_EiJ")    
        
        cclines = paste0("CC", sprintf(1:1000, fmt = "%003d"))
        names(cclines) = cclines
        all.lines = c(names(founders), cclines)
        
        strains = dir(tabledir, full.names=F)
        
        ##    strains = all.lines[strains]
        strains = factor(strains, levels = all.lines)
        strains = sort(strains)
        strains = as.character(strains)
        
        ##    strains = setdiff(strains, "OR559b38V01")
        ##   strains = c("CC001", "CC002")
    
        return(strains)
    }

    inst$getChrs <- function(strain = inst$getStrains()[1])
    {
        straindir = fp(tabledir, strain)
        chrfls    = dir(straindir, pattern = ".*\\.txt\\.tar\\.gz")
        chrs      = gsub(chrfls, pattern = "\\.txt\\.tar\\.gz*", replacement = "")
        
        chrsAll = c(1:19, "X", "Y", "MT")
        chrs = factor(chrs, levels = chrsAll)
        chrs = as.character(sort(chrs))
        return(chrs)
    }

    inst$get.folder <- function(subdir)
    {
        straindir = fp(tabledir, subdir)
        dir.create(straindir, showWarnings = F, recursive = T)
        return(straindir)
    }
    
    inst$getfile <- function(strain, chr, zipped = F)
    {
        straindir = fp(tabledir, strain)
        if(zipped)
        {
            fle = fp(straindir, paste0(chr, ".txt.tar.gz"))
        } else {
            fle = fp(straindir, paste0(chr, ".txt"))
        }
    }
    
    inst$read <- function(strain, chr, zipped=T, select=NULL)
    {
        readargs = list(stringsAsFactors= T, skip=1, sep = inst.sep, showProgress = F)
        fle = inst$getfile(strain, chr, zipped = zipped)
        readargs$select = select ##must be on its own line or else it actually stores a null
        readargs$input = ifelse(zipped, paste0("tar -xOzf ", fle), fle)

        print(paste(strain, chr))
        df = do.call(fread, readargs)

        return(df)
    }

    inst$write <- function(strain, chr, df, zipped = T)
    {
        print(paste0("writing ", chr, " ", strain))
        headerlines = paste0("#strain:", strain, ", chr:",chr,"\n")
        
        straindir = fp(tabledir, strain)
        dir.create(straindir, showWarnings = F, recursive = T)
        fle     = inst$getfile(strain, chr, zipped = F)
        if(!is.null(headerlines))
        {
            headerlines = paste0(headerlines, paste(colnames(df), collapse = inst.sep), "\n")
            cat(headerlines, file=fle)
        }
        fwrite(df, file = fle, row.names = F, sep = inst.sep, append =T, showProgress = F)

        if(zipped)
        {
            fle.zip = inst$getfile(strain, chr, zipped = T)
            
            old.wd = getwd()
            setwd(dirname(fle))

            command = paste("tar -czf ", basename(fle.zip),  basename(fle))
            system(command)
            setwd(old.wd)
            system(paste0("rm ", fle))
        }
        
    }


    return(inst)
}



    ## command = paste0("tar -xzvf ", chrfl.long.tgz, " -C ", straindir)
    ## print(command)
    ## system(command)
    
    ## if(zipped)
    ## {
    
    ## }
