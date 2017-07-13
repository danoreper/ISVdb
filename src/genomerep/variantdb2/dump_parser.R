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


    dblite$.callparse <- function(parseFunc, strain1, strain2=NULL, chr, type)
    {
        df = dblite$read(strain1 = strain1, strain2=strain2, chr=chr, type=type)
        if(is.null(strain2))
        {
            out = parseFunc(strain1, chr, df)
        } else {
            out = parseFunc(strain1, strain2, chr, df)
        }
        return(out)
    }
    
    dblite$read <- function(strain1, strain2=NULL, chr, type, storeKnownFields = T, phased=NULL)
    {
        if(is.null(strain2)&&!is.null(phased))
        {
            stop("phase only has meaning if we are crossing two strains")
        }

        if(!is.null(strain2) &is.null(phased))
        {
            stop("assign a value to 'phased' if crossing two strains.")
        }
        
        if(is.null(strain2))
        {
            inst = dblite[[type]]
            df =  inst$read(strain1, chr, storeKnownFields)
        } else {
            funcname = paste0(".read", type, "Pair")
            df   = dblite[[funcname]](strain1, strain2, chr, storeKnownFields, phased) 
        }
        return(df)
    }

    dblite$.readgenotypePair <- function(strain1, strain2, chr, storeKnownFields, phased = NULL)
    {
        inst = dblite[[paste0("genotype","Sampling")]]
        df1 =  inst$read(strain1, chr, storeKnownFields=F)
        df2 =  inst$read(strain2, chr, storeKnownFields=F, select =c("variant_id", "transcript_name", "allele", "consequence", "prob"))
        setnames(df1, c("allele","consequence"),  c("allele1", "consequence1"))
        setnames(df2, c("allele","consequence"),  c("allele2", "consequence2"))
        
        tojoinback = unique(df1[,c("variant_id", "transcript_name", "pos","gene_name")])
        df1[,pos:=NULL]
        df1[,gene_name:=NULL]
        
        df = df1[df2, on=c("variant_id", "transcript_name"), allow.cartesian=T]
        df[,prob:=prob*i.prob]
        df[,i.prob:=NULL]

        for(cname in c("allele1", "allele2", "consequence1", "consequence2"))
        {
            df[[cname]] = as.character(df[[cname]])
        }

        if(is.null(phased) | !phased)
        {
            correctOrder = df$allele1 <= df$allele2
            
            tmp = df$allele1
            df$allele1[!correctOrder] = df$allele2[!correctOrder]
            df$allele2[!correctOrder] = tmp[!correctOrder]
            
            tmp = as.character(df$consequence1)
            df$consequence1[!correctOrder] = df$consequence2[!correctOrder]
            df$consequence2[!correctOrder] = tmp[!correctOrder]
            
        }

        df = df[,list(consequence1 = consequence1[1], consequence2 = consequence2[1], prob=sum(prob)),
                by=c("variant_id", "transcript_name", "allele1", "allele2")]
        
        df = tojoinback[df, on=c("variant_id","transcript_name"), allow.cartesian=T]
        
        if(storeKnownFields)
        {
            df[,strain1:=strain1]
            df[,strain2:=strain2]
            df[,chr:= chr]
        }
        ordering = c("strain1", "strain2", "variant_id", "chr", "pos", "transcript_name", "gene_name", "allele1", "allele2", "consequence1", "consequence2", "prob")
        newInds = match(ordering, colnames(df))
        newInds = na.omit(newInds)
        df = df[,colnames(df)[newInds],with=F]
        return(df)
    }

    dblite$.readdiplotypePair <- function(strain1, strain2, chr, storeKnownFields, phased = NULL)
    {
        inst = dblite[[paste0("diplotype","Sampling")]]
        df1 =  inst$read(strain1, chr, storeKnownFields=F)
        df2 =  inst$read(strain2, chr, storeKnownFields=F)
        
        setnames(df1, "founder", "founder1")
        setnames(df2, "founder", "founder2")
        
        tojoinback = unique(df1[,c("variant_id", "pos")])
        
        df = df1[df2, on=c("variant_id","gene_name"), allow.cartesian=T]
        df[,prob:=prob*i.prob]
        df[,i.prob:=NULL]
        
        if(is.null(phased) | !phased)
        {
            df[,tmp:=founder1]
            df[,founder1:=pmin(as.character(founder1), as.character(founder2))]
            df[,founder2:=pmax(as.character(tmp), as.character(founder2))]
        }
        df = df[,list(prob=sum(prob)), by=c("variant_id", "gene_name", "founder1", "founder2")]
        df = df[tojoinback, on="variant_id", allow.cartesian=T]
        
        if(storeKnownFields)
        {
            df[,strain1:=strain1]
            df[,strain2:=strain2]
            df[,chr:= chr]
        }
        ordering = c("strain1", "strain2", "variant_id", "chr", "pos", "gene_name", "founder1", "founder2","prob")
        newInds = match(ordering, colnames(df))
        newInds = na.omit(newInds)
        df = df[,colnames(df)[newInds],with=F]
        return(df)
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
    
    inst$read <- function(strain, chr, zipped=T, select=NULL, storeKnownFields = T)
    {
        readargs = list(stringsAsFactors= T, skip=1, sep = inst.sep, showProgress = F)
        fle = inst$getfile(strain, chr, zipped = zipped)
        readargs$select = select ##must be on its own line or else it actually stores a null
        readargs$input = ifelse(zipped, paste0("tar -xOzf ", fle), fle)

        print(paste(strain, chr))
        df = do.call(fread, readargs)
        if(storeKnownFields)
        {
            df$strain = strain
            df$chr     = chr
        }
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
