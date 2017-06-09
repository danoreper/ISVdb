##Highest level functions for building variant dbs. 
library(IRanges)
library(data.table)
library(reshape2)
library("Biostrings")

source("./loadParams.R")
source("./genomerep/cc_founderProbs.R")
source("./genomerep/filtervcf.R")
source("./genomerep/buildGenomeData2.R")
source("./stringutils.R")
source("./genomerep/cc_founderProbs.R")
source("./genomerep/variantdb2/dump_parser.R")


buildVariantDb = new.env(hash=T)

getIterAccumulator <- function(batchSize=1, mem.gb = 20)
{
    if(prop$onCluster)
    {
        bsubCommand = bsub$get.default.killdevil.bsub(numProcessPerNode = 1, memoryLimit.GB = mem.gb, queue="day")
        accum = bsub$get.bsub.accumulator("./genomerep/variantdb2/buildVariantDb.R", bsubCommand, batchSize=batchSize)
    } else {
        
        corez = 1 #prop$mnp$mc.cores
        batchSize = corez*100
        accum = bsub$get.mc.accumulator(mc.cores= corez)
    }
    return(accum)
}


##primary entry point:
##build all the variant DBs, exons for now, imprinted and whole genome later
buildVariantDb$buildDbs <- function(genomeData=buildGenomeData$buildAllData()) 
{
    .extend <- function(x, upstream=0, downstream=0)     
    {
        if (any(strand(x) == "*"))
        {
            warning("'*' ranges were treated as '+'")
        }
        on_plus <- strand(x) == "+" | strand(x) == "*"
        new_start <- start(x) - ifelse(on_plus, upstream, downstream)
        new_end <- end(x) + ifelse(on_plus, downstream, upstream)
        ranges(x) <- IRanges(new_start, new_end)
        trim(x)
    }
    
    karyotype = genomeData$karyotype
##    for(dbtype in c("full", "exon")) ##full, ##, "imprinted")) TODO full and imprinted once disk is available.
    for(dbtype in c("exon", "full")) ##full, ##, "imprinted")) TODO full and imprinted once disk is available.
    {
        print(paste0("Building dbtype:", dbtype))
        ranges = NULL
        if(dbtype == "exon")
        {
            ##TODO check to see if this is already in genome data
            print("building exons")
            ranges = buildGenomeData$getExons()
            print("done building exons")
        } 

        ##pad relevant ranges by some amount if not using full db.
        if(!is.null(ranges))
        {
            thepadding =  prop$variantdb[[dbtype]][["padding"]]
            ranges = .extend(ranges, upstream = thepadding, downstream = thepadding)
        }
                                             
        dbDir   = dat(fp(prop$variantdb$dbfolder, prop$variantdb[[dbtype]]$name))
        db = db_builder$get.db.lite(dbDir)
        
        
        buildVariantDb$buildSingleVariantDb(db = db,
                                            ranges = ranges,
                                            transcriptsFile = NULL,
                                            rebuildVCF      = prop$variantdb[[dbtype]]$rebuildVCF,
                                            rebuildFounder  = prop$variantdb[[dbtype]]$rebuildFounder,
                                            rebuildCC       = prop$variantdb[[dbtype]]$rebuildCC,
                                            limit           = prop$variantdb$var_limit)
    }
}

buildVariantDb$.createBedForChr <- function(db, karyotype, ranges, chr)
{
    len = karyotype$len[karyotype$chrname == chr]
    gr = c(buildGenomeData$buildGr(start=1, end = len, chr = chr, strand = "+"),
           buildGenomeData$buildGr(start=1, end = len, chr = chr, strand = "-"))
    
    if(!is.null(ranges))
    {
        gr = intersect(gr, ranges)
    }
    bed.dir = fp(db$get.temp.dir("bed"))
    chrbed = fp(bed.dir, paste0("chr_",chr,".bed"))
    buildGenomeData$writeToBed(gr, chrbed)
    return(chrbed)
}

##build a single variant db.
##dbname: name for the db
##chr: the chromosome for the db.
##karyotype: a data frame providing chromosome lengths per chromosome.
##filter.bed: a bed file restricting the regions of interest.
##initialVariantId: the startId to build a variantDB. This is typically used
## when we loop over chromosomes, and want a unique id for every variant, so start with
## the last chromosome built max variant id + 1.
## transcriptsFile: a list of transcripts of interest. TODO delete me.
## pythonlogfile: location for the python log for the vcf parser.
## rebuildVCF: boolean indicating whether the vcf for this db has already been built. If it has, dont bother rebuilding. Primaruily for testing.
## rebuildFounder: boolean indicating whether the VCF file has been parsed and a founder database for it already built.
## rebuildCC: boolean indicating whether the CC database part has already been built.
## limit: the maximum number of variants to store in the db. PRimarily for testing.

buildVariantDb$buildSingleVariantDb <- function(db,
                                                chrsToBuild     = prop$variantdb$chr_range,
                                                founders        = read.table(dat(prop$genome$foundersMap),header=T, sep=",")$founder,
                                                vcfdir           = dat(prop$genome$vcfDir),
                                                rilHaplotypeProbsDir = dat(prop$genome$rilHaplotypeProbsDir),
                                                karyotype       = buildGenomeData$getKaryotype(referenceFile = dat(prop$genome$karyotype)),
                                                ranges          = NULL,
                                                transcriptsFile = NULL,
                                                limit           = prop$variantdb$var_limit,
                                                rebuildVCF      = prop$variantdb[[dbtype]]$rebuildVCF,
                                                rebuildFounder  = prop$variantdb[[dbtype]]$rebuildFounder,
                                                rebuildCC       = prop$variantdb[[dbtype]]$rebuildCC)                                                
{
    if(is.na(chrsToBuild))
    {
        chrsToBuild = karyotype$chrname
    }
    ##build an accumulator for parallelizing jobs
    accum = getIterAccumulator()
    ##accum = bsub$get.stub.accum()

    ##Rebuild founder genotypes
    pracma::tic()
    if(rebuildFounder)
    {
        print("rebuilding founder data")
        accum$init(func = buildVariantDb$buildFounderGenotypes,
                   otherGlobals = list(db              = db,
                                       founders        = founders,
                                       vcfdir          = vcfdir,
                                       karyotype       = karyotype,
                                       ranges          = ranges,
                                       transcriptsFile = transcriptsFile,
                                       limit           = limit,
                                       rebuildVCF      = rebuildVCF))
        for(chr in chrsToBuild)
        {
            accum$addCall(list(chr = chr))
        }
        accum$runAll()
        print("updating variant_id to be unique across chromosomes")
        buildVariantDb$updateFounderVariantIDs(db, founders)
    }
    pracma::toc(); print("done parsing vcf")

    
    ##Rebuild CC info
    ##accum = bsub$get.stub.accum()
    ##chrsToBuild = "MT"
    accum = getIterAccumulator(mem.gb=25)##getIterAccumulator(mem.gb = 20, batchSize = 10)
    pracma::tic(); print("began building ccdb")
    if(rebuildCC)
    {
        accum$init(func = buildVariantDb$build_CC_info,
                   otherGlobals = list(db = db,
                                       founders = founders,
                                       rilHaplotypeProbsDir = rilHaplotypeProbsDir,
                                       karyotype = karyotype))
        
        for(i in 1:length(chrsToBuild))
        {
            chr = chrsToBuild[[i]]
            accum$addCall(list(chr = chr))
        }
        accum$runAll()
        print("finished building ccdf")
    }
}

buildVariantDb$updateFounderVariantIDs <- function(db, founders)
{

    chrsToCheck =  db$genotype$getChrs()
    maxv = unlist(db$iterate("genotype", parseFunc = function(strain,chr, df){max(df$variant_id)},
                             strain1s = "C57BL6J", chrs = chrsToCheck, select = "variant_id"))
    
    offset = c(0, cumsum(maxv)[-length(maxv)]) ##the first element doesnt need an offset.. shift all of them over
    names(offset) =chrsToCheck
    
    helper.func <- function(strain, chr,df)
    {
        df[,variant_id := variant_id + offset[chr]]
        db$genotype$write(strain,chr,df)
    }

    accum = getIterAccumulator(batchSize=8, mem.gb = 6) ##bsub$get.stub.accum() 
    db$iterate("genotype", parseFunc = helper.func, accum = accum, strain1s = founders)
}


##TODO: 
assignMax <- function(df, bycolz)
{
    df[,isMax:=NULL]
    df[,isMax := as.integer(prob==max(prob)), bycolz]
    allcolz = colnames(df)
    ##move ismax to location just after prob 
    allcolz = setdiff(colnames(df), "isMax")
    allcolz = util$insertAtIndex("isMax", which(allcolz=="prob")+1, allcolz) 
    setcolorder(df, allcolz)

    ## pracma::tic()
    ## df[,isMax := as.integer(1:.N==which.max(prob)), by=c("variant_id", "transcript_name")]
    ## pracma::toc()
    
    ## pracma::tic()
    ## df[,isMax3 := as.integer(frank(-prob, ties.method = "min")==1), by=c("variant_id", "transcript_name")]
    ## pracma::toc()
}


##TODO: pass in a founder or a founder file?
##Wrap up R data into a python command that calls the vcf_parser
buildVariantDb$buildFounderGenotypes <- function(db,
                                                 chr,
                                                 founders,
                                                 vcfdir,
                                                 karyotype,
                                                 ranges = NULL,
                                                 transcriptsFile = NULL,
                                                 limit = NA,
                                                 rebuildVCF = T)
{
    ##Filter vcf to only include chromosome. If ranges are specified, filter the vcf to only those ranges.
    print("filtering vcf")

    vcftargetdir   = db$get.temp.dir(fp("vcf",paste0("chr_",chr)))
    pracma::tic()
    if(rebuildVCF)
    {
        ##the default args to filtervcf which change depending on whether a set of exons is passed in.
        filter.args = list(vcfDirIn       = vcfdir,
                           vcfDirOut      = vcftargetdir,
                           tmpdir         = prop$tmpdir, 
                           vcfExtensionIn =  ".vcf.gz", 
                           founders       = founders, 
                           mc.cores       = 1)

        if(!is.null(ranges))
        {
            chrbed = buildVariantDb$.createBedForChr(db, karyotype, ranges, chr)
            filter.args[["bedfile"]] = chrbed
        } else {
            filter.args[["chr"]] = chr
        }
        do.call(filterVCF$filterVcfsInDir, filter.args)
    } 
    pracma::toc() 
    print("done filtering vcf")
    

    ##Given a vcf filtered to only include the relevant chromosomes, parse it and build tables of founder genotypes.
    pracma::tic()
    print("parsing vcf")
    vcfParser  = "nohup python -m genomerep.variantdb2.vcf_parser "

    
    transcriptString = ""
    if(!is.null(transcriptsFile))
    {
        transcriptString     = paste0("-t ", transcriptsFile)
    }
    
    limitString = ""
    if(!is.na(limit))
    {
        limitString = paste0("-l ", limit)
    }

    pythonlogfile = fp(db$get.temp.dir("vcfparser_logs"), paste0("chr_", chr, ".log"))
    pythonlogstring = paste0(" 2> ", pythonlogfile)

    foundersFile = fp(vcftargetdir, "founder.csv")
    founders = data.table(founders)
    fwrite(data.table(founders), foundersFile, col.names=F, row.names=F, quote=F)
    founder.genotypes.file = fp(db$get.temp.dir("founders"), paste0("chr_", chr, ".txt"))
    command = paste(vcfParser,
                    vcftargetdir,
                    foundersFile,
                    founder.genotypes.file,
                    transcriptString,
                    limitString,
                    pythonlogstring)
    
    print(command)
    system(command)
    pracma::toc()
    print("done parsing vcf")

    ## Split up the big founder genotypes file, and write them each out to their own file.
    founder.genotypes = fread(founder.genotypes.file, stringsAsFactors = T)
    assignMax(founder.genotypes, bycolz = c("strain","variant_id","transcript_name"))

    print("about to split files")
    ##do this one at a time to avoid blowing up memory
    for(astrain in levels(founder.genotypes$strain))
    {
        print(paste0("splitting ", astrain))
        adf = founder.genotypes[strain == astrain]
        adf[,chr:=NULL]
        adf[,strain:=NULL]
        assignMax(adf, bycolz = c("variant_id","transcript_name"))
        db$genotype$write(astrain, chr, adf) 
    }

    ##remove the redundant version, avoid taking up space.
    ##unlink(founder.genotypes.file)
    
##    curdir = dirname(founder.genotypes.file)      
##    split.perfile = paste0("awk -F, \'{print >> (\"",  curdir, "/\"$3\".txt\")}\' ",founder.genotypes.file)


}






buildVariantDb$build_CC_info <- function(db,
                                         chr,
                                         founders,
                                         rilHaplotypeProbsDir,
                                         karyotype)
                       
{
    
    getDiplotypeTable <- function(probFile, variant_id, chr, pos, genez.df)
    {
        probData   = buildFounderProbs$build(probFile, founders,karyotype, chr=chr)
        dp = buildFounderProbs$hapProbPerPosition(variant_id, chr, pos, probData)
        
        ##build mapping from a diplotype string to the componenet haplotypes;
        ##regex operation can be surprisingly slow, so we are caching the component haplotypes
        ##once, and then using that lookup for every variant.
        alldiplo = colnames(dp)
        id1 = sub("\\..*","",alldiplo)
        id2 = sub(".*\\.","",alldiplo )
        map1 = id1
        names(map1) = alldiplo
        map2 = id2
        names(map2) = alldiplo
        
        dp = convert_GT_to_DB_format(gt = dp, genez.df = genez.df, map1, map2)
        return(dp)
    }

    convert_GT_to_DB_format <- function(gt, genez.df, map1, map2) 
    {
        melted = data.table::melt(gt,id=c("variant_id"), variable.name = "diplotype", value.name="prob")
        melted = melted[melted$prob>=prop$variantdb$truncate_cc_prob,]
        rm(gt)
        
        chars = as.character(melted$diplotype)
        id1 = map1[chars]
        id2 = map2[chars]
        
        melted$founder_1 = pmin(id1, id2)
        melted$founder_2 = pmax(id1, id2)
        melted = melted[,c("variant_id", "founder_1", "founder_2", "prob"), with = F]
        setkey(melted, "variant_id")
        melted = genez.df[melted, allow.cartesian = T]
        setcolorder(melted, c("variant_id", "pos", "founder_1", "founder_2", "prob", "gene_name"))
        assignMax(melted, c("variant_id", "gene_name"))
        return(melted)
    }


    getFounderSamplingTable <- function(db, chr, founders)
    {
        allfounders = list()
        for(founder in founders)
        {
            df = db$genotype$read(strain = founder, chr = chr)
            df[,strain:=founder]
            allfounders = util$appendToList(allfounders, df)
        }
        allfounders = rbindlist(allfounders)
        ## ##TODO: build genotype sampling files, use them.
        ## allfounders.file = fp(db$get.temp.dir(fp("founders")), paste0("chr_", chr,".txt"))
        ## allfounders = fread(allfounders.file)

        allfounders[,isMax:=NULL]
        gt.founder.sampling = getSamplingTable(allfounders, c("allele", "consequence"), c("strain","variant_id","transcript_name","allele"))
        return(gt.founder.sampling)
        ##return(dfs)
    }
    
    getGenotypeTable <- function(dp, gt.founder.sampling)
    {
        dp[,pos:=NULL]
        dp[,gene_name:=NULL]
        dp = unique(dp)

        setkey(dp, "variant_id", "founder_1")
        setkey(gt.founder.sampling, "variant_id", "strain")
        gt1 = gt.founder.sampling[dp]
        rm(dp)
        gc()
        
        gt1[,prob:=prob*i.prob]
        gt1[,i.prob:=NULL]
        gt1[,pos:=NULL]
        gt1[,gene_name:=NULL]
        gt1[,strain:=NULL]
        
        setkey(gt1, "variant_id", "founder_2", "transcript_name")
        setkey(gt.founder.sampling, "variant_id", "strain", "transcript_name")
        gt2 = gt.founder.sampling[gt1]
        
        gt2[,prob:=prob*i.prob]
        gt2[,i.prob:=NULL]
        gt2[,strain:=NULL]
        
        setnames(gt2,
                 old = c("allele", "i.allele", "consequence", "i.consequence"),
                 new = c("allele_1", "allele_2", "consequence_1", "consequence_2"))
        
        setcolorder(gt2, c("variant_id", "pos", "allele_1", "allele_2", "prob", "gene_name", "transcript_name", "consequence_1", "consequence_2"))
        assignMax(gt2, bycolz=c("variant_id", "transcript_name"))
        return(gt2)
    }

    
    getSamplingTable <- function(allfounders, meltover, sumover)
    {
        colz1 = setdiff(colnames(allfounders), paste0(meltover, "_2"))
        colz2 = setdiff(colnames(allfounders), paste0(meltover, "_1"))
        df1 = allfounders[,colz1,with=F]
        df2 = allfounders[,colz2,with=F]
        
        setnames(df1, old = paste0(meltover, "_1"), new = meltover)
        setnames(df2, old = paste0(meltover, "_2"), new = meltover)
        
        x = rbind(df1, df2)
        y = x[,c(head(.SD,1), list(prob=sum(prob))), by = sumover, .SDcols = setdiff(colnames(x), c(sumover, "prob"))]

        y[,prob := prob*.5]
        return(y)
    }          


    
    ##Get founder genotype sampling table
    gt.founder.sampling = getFounderSamplingTable(db, chr=chr, founders = founders)
    
    ##unique genes
    genez.df = unique(db$genotype$read("C57BL6J", chr, zipped = T, select = c("variant_id", "pos","gene_name")))
    setkey(genez.df, "variant_id")
    
    ##unique variants
    var.df     = unique(genez.df[,c("variant_id", "pos")])

    
    probsFiles = stringutils$getFilesWithExtension(rilHaplotypeProbsDir, ".csv")
    cc_limit   = ifelse(is.na(prop$variantdb$cc_limit), length(probsFiles), prop$variantdb$cc_limit)
    for(i in 1:cc_limit)
    {
        pracma::tic()
        probFile = probsFiles[i]
        print(probFile)
        
        strainNameAndDescent = buildFounderProbs$parseProbFilename(probFile)
        strainName = strainNameAndDescent$strain[[1]]

        ##get the diplotype table, write it to file
        print("getting diplotypes")
        dp = getDiplotypeTable(probFile, var.df$variant_id, chr, var.df$pos, genez.df)
        db$diplotype$write(strain = strainName, chr = chr, df = dp, zipped = T)

        
        ##get the diplotype sample table, write it to file
        print("getting diplo sampling") 
        dp[,isMax:=NULL]
        dp.sampling  = getSamplingTable(dp, meltover = "founder", sumover = c("variant_id", "founder","gene_name"))
        db$diplotypeSampling$write(strain = strainName, chr = chr, df = dp.sampling, zipped = T)
        rm(dp.sampling)
        gc()
        
        ##get the genotype table, write it to file
        print("get genotypes")
        gt = getGenotypeTable(dp, gt.founder.sampling)
        db$genotype$write(strain = strainName, chr = chr, df = gt, zipped = T)
        rm(dp)
        gc()
        
        ##get the genotype sampling table, write it to file
        print("get genotype sampling")
        gt[,isMax:=NULL]
        gt.sampling = getSamplingTable(gt, meltover = c("allele","consequence"), c("variant_id","transcript_name","allele"))
        db$genotypeSampling$write(strain = strainName, chr = chr, df = gt.sampling, zipped = T)
        
        rm(gt)
        rm(gt.sampling)
        gc()
        pracma::toc()
        
    }
    
    ##add founder diplotypes-- trivial, but necessary for backcrosses
    for(founder in founders)
    {
        print(founder)

        ##write founder diplotypes to file
        dp.founder = copy(genez.df)
        dp.founder$founder_1 = founder
        dp.founder$founder_2 = founder
        dp.founder$prob      = 1
        dp.founder$isMax     = 1
        setcolorder(dp.founder, c("variant_id", "pos", "founder_1", "founder_2", "prob", "isMax", "gene_name"))
        db$diplotype$write(strain=founder, chr=chr, df = dp.founder, zipped = T)

        ##write founder sampling diplotypes to file
        dp.founder[,isMax:=NULL]
        dp.founder.sampling  = getSamplingTable(dp.founder, meltover = "founder", sumover = c("variant_id", "founder","gene_name"))
        db$diplotypeSampling$write(strain = founder, chr = chr, df = dp.founder.sampling, zipped = T)
    }

    ##add founder genotype sampling. Not quite as trivial given residual heterozygosity, but we already computed it above
    for(founder in founders)
    {
        gt.sampling = gt.founder.sampling[strain == founder]
        gt.sampling[,chr:=NULL]
        gt.sampling[,strain:=NULL]
        db$genotypeSampling$write(strain = founder, chr = chr, df = gt.sampling, zipped = T)
    }
    
}





