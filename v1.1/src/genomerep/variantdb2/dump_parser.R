source("./bsub.R")
source("./loadParams.R")

table_manip$get.db.lite(dbdir)
{
    dblite = new.env(hash=T)
    dblite$genotype  = table_manip$getInstance(fp(dbdir, "genotype"))
    dblite$diplotype = table_manip$getInstance(fp(dbdir, "diplotype"))
    dblite$temp      = fp(dbdir, temp)##table_manip$getInstance(fp(dbdir, "temp") #fp(dbdir, temp)
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
    
    return(dblite)
}

table_manip = new.env(hash=T)
##TODO: property?
##table_manip$strainsdir = dat(fp("dump_100bp_2017-03-11/genotype"))

table_manip$getInstance <- function(tabledir)
{
    dir.create(tabledir, showWarnings = F, recursive = T)
    force(tabledir)
    inst = new.env(hash = T)

    inst$clear <- function()
    {
        
    }
    
    inst$iterate <- function(parseFunc,
                             strains = inst$getStrains(),
                             chrs = inst$getChrs(strain = strains[1]),
                             zipped = T,
                             accum = bsub$get.mc.accumulator(mc.cores = 1),
                             iterator = F)    
    {
        outs = list()
        force(parseFunc)

        accum$init(func = inst$.callparse, otherGlobals= list(parseFunc = parseFunc, zipped = zipped))
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
        fle = inst$getfile(strain, chr, zipped = zipped)        
        if(zipped)
        {
            command = paste0("tar -xOzf ", fle)
            df = fread(input=command, stringsAsFactors=T, skip =1, select = select)
        } else { 
            df = fread(fle, stringsAsFactors = T, skip = 1, select = select)
        }
        
        ##TODO remove
        ##  col.names = c("variant_id", "pos", "allele_1", "allele_2", "prob", "is_max", "consequence_1", "consequence_2", "gene_name", "transcript_name")
        ## setnames(df, old = colnames(df), new=col.names)
    }

    inst$write <- function(strain, chr, df, zipped = T, headerlines=NULL)
    {
        straindir = fp(table_manip$strainsdir, strain)
        dir.create(straindir, showWarnings = F, recursive = T)
        fle     = inst$getfile(strain, chr, zipped = F)
        if(!is.null(headerlines))
        {
            cat(headerlines, file=fle)
        }
        fwrite(df, file = fle, row.names = F, sep = "\t", append =T)

        if(zipped)
        {
            fle.zip = inst$getfile(strain, chr, zipped = T)
            
            old.wd = getwd()
            setwd(dirname(fle))

            command = paste("tar -czvf ", basename(fle.zip),  basename(fle))
            system(command)
            setwd(old.wd)
            system(paste0("rm ", fle))
        }
        
    }


    inst$.callparse <- function(parseFunc, strain, chr, zipped = T)
    {
        df =  inst$read(strain, chr, zipped = zipped)
        out = parseFunc(strain, chr, df)
        return(out)
    }

    return(inst)
}



    ## command = paste0("tar -xzvf ", chrfl.long.tgz, " -C ", straindir)
    ## print(command)
    ## system(command)
    
    ## if(zipped)
    ## {
    
    ## }
