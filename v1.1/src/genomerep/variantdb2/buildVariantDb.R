##Highest level functions for building variant dbs. 
library(RMySQL)
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
source("./genomerep/variantdb/dump_parser.R")



buildVariantDb = new.env(hash=T)

##primary entry point:
##build all the variant DBs, exons for now, imprinted and whole genome later
buildVariantDb$buildDbs <- function(genomeData=buildGenomeData$buildAllData) 
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
    for(dbtype in c("exon")) ##full, ##, "imprinted")) TODO full and imprinted once disk is available.
    {
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
                                             
        dbDir   = fp(prop$variantdb$dbfolder, prop$variantdb[[dbtype]]$name)
        db = table_manip$get.db.lite(dbDir)
        
        
        buildVariantDb$buildSingleVariantDb(db = db,
                                            ranges = ranges,
                                            transcriptsFile = NULL,
                                            rebuildVCF      = prop$variantdb[[dbtype]]$rebuildVCF,
                                            rebuildFounder  = prop$variantdb[[dbtype]]$rebuildFounder,
                                            rebuildCC       = prop$variantdb[[dbtype]]$rebuildCC,
                                            limit           = prop$variantdb$var_limit)
    }
}

buildVariantDb$.createBedForChr(db, karyotype, ranges, chr)
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
                                                rilHaplotypeProbsDir = prop$genome$rilHaplotypeProbsDir,
                                                karyotype       = NULL,
                                                ranges          = NULL,
                                                transcriptsFile = NULL,
                                                limit           = prop$variantdb$var_limit,
                                                rebuildVCF      = prop$variantdb[[dbtype]]$rebuildVCF,
                                                rebuildFounder  = prop$variantdb[[dbtype]]$rebuildFounder,
                                                rebuildCC       = prop$variantdb[[dbtype]]$rebuildCC)
                                                
{

    if(is.null(karyotype))
    {
        karyotype = buildGenomeData$getKaryotype(referenceFile = dat(prop$genome$dnaReferenceFile))
    } 

    if(is.na(chrsToBuild))
    {
        chrsToBuild = karyotype$chrname
    }
    
    for(chr in chrsToBuild)
    {
        print(paste0("building chr: ",chr))

        pracma::tic()
        if(rebuildFounder)
        {
            buildVariantDb$buildFounderGenotypes(db              = db,
                                                 chr             = chr,
                                                 founders        = founders
                                                 vcfdir          = vcfdir,
                                                 karyotype       = karyotype,
                                                 ranges          = ranges,
                                                 transcriptsFile = transcriptsFile,
                                                 limit           = limit,
                                                 rebuildVCF      = rebuildVCF)
        }
        pracma::toc(); print("done parsing vcf_parser")
        
        pracma::tic(); print("began building ccdb")
        if(rebuildCC)
        {
            ##TODO read in only B6 rather than taking unique here
            var.df = unique(db$genotype$read("founders", chr, zipped = F, select = c("variant_id", "pos")))

            buildVariantDb$build_CC_diplotypes(db,
                                               chrom = chr,
                                               founders = founders,
                                               rilHaplotypeProbsDir = rilHaplotypeProbsDir,
                                               var.df$variant_id,
                                               pos = inp$pos)
            
            print("finished building ccdf")
        }
    }
    ##TODO: do this after all dbs are built.
    ##initialVariantId = buildVariantDb$.getInitVariantId(db.manip, karyotype, chrsToBuild)

}

##TODO: pass in a founder or a founder file?
##Wrap up R data into a python command that calls the vcf_parser
buildVariantDb$buildFounderGenotypes(db,
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
                           ##        tmpdir         = tmpdir, 
                           vcfExtensionIn =  ".vcf.gz", 
                           founders       = founders, 
                           mc.cores       = 1)

        if(!is.null(ranges))
        {
            chrbed = buildVariantDb$.createBedForChr(db, karyotype, ranges, chr)
            filter.args[["bedfile"]] = chrbed
            do.call(filterVCF$filterVcfsInDir, filter.args)
        } else {
            filter.args[["chr"]] = chr
        }
    } 
    pracma::toc() 
    print("done filtering vcf")
    

    ##Given a vcf filtered to only include the relevant chromosomes, parse it and build tables of founder genotypes.
    pracma::tic()
    print("parsing vcf")
    vcfParser  = "nohup python -m genomerep.variantdb.vcf_parser "

    
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
    fwrite(founders, foundersFile, col.names=F, row.names=F, quote=F)

    founder.genotypes.file = fp(db$get.temp.dir("founders"), paste0("chr_", chr, ".log"))
    command = paste(vcfParser,
                    vcftargetdir,
                    foundersFile,
                    founder.genotypes.file,
                    transcriptString,
                    limitString,
                    chrstring,
                    pythonlogstring)
    
    pracma::toc()
    print("done parsing vcf")
    
    print(command)
    system(command)
    Sys.sleep(5)
    return(founder.genotypes.file)
}

buildVariantDb$build_CC_diplotypes <- function(db,
                                               chrom,
                                               founders,
                                               rilHaplotypeProbsDir,
                                               variant_id,
                                               pos,
                                               karyotype)
{
    probsFiles = stringutils$getFilesWithExtension(rilHaplotypeProbsDir, ".csv")
    
    cc_limit = prop$variantdb$cc_limit
    if(is.na(cc_limit))
    {
        cc_limit = length(probsFiles)
    }



    ##strainName = paste(strainNameAndDescent$strain, strainNameAndDescent$descent, sep="_")
    for(probFile in probsFiles[1:cc_limit])
    {
        strainNameAndDescent = buildFounderProbs$parseProbFilename(probFile)
        strainName = strainNameAndDescent$strain
        probData   = buildFounderProbs$build(probFile, founders,karyotype)
        gt = buildFounderProbs$hapProbPerPosition(variant_id, chrom, pos, probData)
        db_version_of_gt = buildVariantDb$convert_GT_to_DB_format(gt = gt)

        ##TODO: zip?
        db$diplotype$write(strain = strainName, chr = chrom, df = db_version_of_gt, zipped = T,
                           headerlines = paste0("strain:", strainName, ", chr:",chrom))
    }
}


buildVariantDb$convert_GT_to_DB_format <- function(gt) 
{
    melted = data.table::melt(gt,id=c("variant_id"), variable.name = "diplotype", value.name="prob")
    melted = melted[melted$prob>=prop$variantdb$truncate_cc_prob,]
    rm(gt)
    chars = as.character(melted$diplotype)
    melted$founder_1 = pmin(id1, id2)
    melted$founder_2 = pmax(id1, id2)
    melted = melted[,c("variant_id", "founder_1", "founder_2", "prob"), with = F]
    return(melted)
}




    
##TODO reimplement
## ##gets an initial unused variant id for first chromosome in chrsTo build by iterating over all chromosomes NOT to rebbuild
## buildVariantDb$.getInitVariantId <- function(db.manip, karyotype, chrsToBuild)
## {
##     print(chrsToBuild)
##     alreadyBuiltChrs = setdiff(karyotype$chrname, chrsToBuild)
##     if(length(alreadyBuiltChrs)==0)
##     {
##         return(0)
##     }

##     func = function(strain, chr, df)
##     {
##         return(max(df$variant_id))
##     }
    
##     maxs = db.manip$genotype$iterate(func,
##                                       strains = "founders",
##                                       chrs = alreadBuiltChrs,
##                                       zipped = T)
##     initID = max(maxs)+1
##     return(initID)
## }

##  jobBounds = util$getIndexGroupsForLen(nrow(variantFrame), 1000000)
