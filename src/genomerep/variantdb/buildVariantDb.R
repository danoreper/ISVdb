##Highest level functions for building variant dbs. 

source("./loadParams.R")
source("./genomerep/cc_founderProbs.R")
source("./genomerep/variantdb/add_cc_toVariantDb.R")
source("./genomerep/filtervcf.R")
source("./genomerep/buildGenomeData2.R")

vcfParser  = "nohup python -m genomerep.variantdb.vcf_parser "

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
                                             
        chrsToBuild = prop$variantdb$chr_range
        if(is.na(chrsToBuild))
        {
            chrsToBuild = karyotype$chrname
        }
        
        ##iterate over all chromosomes in karyotype that are not going to be rebuilt to get the max variant id; we want every chromosome that we subsequently build to have non-overlapping variant ids
        initialVariantId = buildVariantDb$.getInitVariantId(karyotype,chrsToBuild, dbtype)

        for(chr in chrsToBuild)
        {
            
            print(paste0("building db: ",dbtype, "_",chr))
            len = karyotype$len[karyotype$chrname == chr]
            
            gr = c(buildGenomeData$buildGr(start=1, end = len, chr = chr, strand = "+"),
                   buildGenomeData$buildGr(start=1, end = len, chr = chr, strand = "-"))

            if(!is.null(ranges))
            {
                gr = intersect(gr, ranges)
            }
            chrbed = fp(prop$tmpdir, paste0(dbtype, "_chr_",chr,".bed"))
            buildGenomeData$writeToBed(gr, chrbed)
            
            dbname = paste0(prop$variantdb[[dbtype]]$name, "_", chr)
            buildVariantDb$buildSingleVariantDb(dbname = dbname,
                                                filter.bed = chrbed,
                                                initialVariantId=initialVariantId,
                                                transcriptsFile = NULL,
                                                pythonlogfile   = fp(paste0("./", dbtype,".log")),
                                                rebuildVCF      = prop$variantdb[[dbtype]]$rebuildVCF,
                                                rebuildFounder  = prop$variantdb[[dbtype]]$rebuildFounder,
                                                rebuildCC       = prop$variantdb[[dbtype]]$rebuildCC,
                                                limit           = prop$variantdb$var_limit)
            
            ##dirty hack to avoid checking if db already exists, but whatever
            initialVariantId = try(buildVariantDb$.getMaxVariantId(dbname) + 1)
            if(class(initialVariantId)=="try-error")
            {
                initialVariantId = NA
            }
        }
    }
}

##build a single variant db.
##dbname: name for the db
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


buildVariantDb$buildSingleVariantDb <- function(dbname,
                                                karyotype  = NULL,
                                                filter.bed = NULL,
                                                initialVariantId = NULL,
                                                transcriptsFile = NULL,
                                                pythonlogfile = NULL,
                                                rebuildVCF = T,
                                                rebuildFounder = T,
                                                rebuildCC = T,
                                                limit = NA)
{
    mc.cores = 1
    
    founders = read.table(dat(prop$genome$foundersMap),header=T, sep=",")$founder
    
    tmpdir          = prop$tmpdir
    dbhost          = prop$variantdb$host
    dbuser          = prop$variantdb$user
    dbpassword      = prop$variantdb$password
    foundersFile    = tempfile("founder", tmpdir, ".csv")
    write.table(founders, foundersFile, col.names=F, row.names=F, quote=F)
    
    vcfExtensionIn = ".vcf.gz"
    dir.create(dat(fp(prop$genome$vcfDir, "filt",dbname)), recursive = T, showWarnings = F)
    vcftargetdir = dat(fp(prop$genome$vcfDir, "filt", dbname))
    print("processing vcf")
    pracma::tic()
    if(rebuildVCF)
    {
        filterVCF$filterVcfsInDir(
          vcfDirIn       = dat(prop$genome$vcfDir),
          vcfDirOut      = vcftargetdir,
          tmpdir         = tmpdir, 
          vcfExtensionIn = vcfExtensionIn, 
          founders       = founders, 
          mc.cores       = mc.cores,
          bedfile        = filter.bed)
    }
    pracma::toc()
    
    
    transcriptString = ""
    if(!is.null(transcriptsFile))
    {
        transcriptString     = paste0("-t ", transcriptsFile)
    }

    initialVariantString = ""
    if(!is.null(initialVariantId))
    {
        initialVariantString = paste0("-i ", initialVariantId)
    }

    limitString = ""
    if(!is.na(limit))
    {
        limitString = paste0("-l ", limit)
    }

    engineString = paste0("-e ", engine)
    
    if(is.null(pythonlogfile))
    {
        pythonlogfile = "./python.log"
    }

    print("parsing vcf")
    pracma::tic()
    pythonlogstring = paste0(" 2> ", pythonlogfile)
    if(rebuildFounder)
    {
        command = paste(vcfParser,
                        vcftargetdir,
                        foundersFile,
                        dbname,
                        dbhost,
                        dbuser,
                        dbpassword,
                        initialVariantString,
                        transcriptString,
                        engineString,
                        limitString,
                        pythonlogstring)
        print(command)
        system(command)
        Sys.sleep(5)
    }
    pracma::toc(); print("done building founder db")

    pracma::tic(); print("began building ccdb")
    if(rebuildCC)
    {
        if(is.null(karyotype))
        {
            karyotype = buildGenomeData$getKaryotype(referenceFile = dat(prop$genome$dnaReferenceFile))
        } 

        rebuild_ccdb(dbname,
                     dat(prop$genome$rilHaplotypeProbsDir),
                     founders,
                     karyotype,
                     dbhost,
                     dbuser,
                     dbpassword)
        print("finished building ccdf")
    }
}


##gets an initial unused variant id for first chromosome in chrsTo build by iterating over all chromosomes NOT to rebbuild
buildVariantDb$.getInitVariantId <- function(karyotype, chrsToBuild,dbtype)
{
    print(chrsToBuild)
    maxVariantId = 0
    for(chr in karyotype$chrname)
    {
        if(! chr %in% chrsToBuild)
        {
            dbname = paste0(prop$variantdb[[dbtype]]$name, "_", chr)
            candMax = try(buildVariantDb$.getMaxVariantId(dbname))
            if(class(candMax) == "try-error")
            {
                next
            } else {
                maxVariantId = max(candMax, maxVariantId)
            }
        }
    }

    initialVariantId = maxVariantId + 1
    return(initialVariantId)
}

buildVariantDb$.getMaxVariantId <- function(db,
                                            dbhost = prop$variantdb$host,
                                            dbuser = prop$variantdb$user,
                                            dbpassword = prop$variantdb$password)
{
    con = dbConnect(MySQL(), user=dbuser, password =dbpassword, dbname=db, host=dbhost)   
    maxid = dbGetQuery(con, "select max(variant_id) from variant;")[1,1]
    dbDisconnect(con)
    return(maxid)
}
