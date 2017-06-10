## Wrapper around the VCFtools command line calls. Used to filter vcfs by strain and genomic region;
## the vcfs we use start out as whole genome and include something like 20 strains, of which
## we only are interested in 7.
## Author: doreper
###############################################################################
source("./stringutils.R")
filterVCF = new.env(hash=T)


##filter a single vcf file to just the founders of interest, in only the region specified by bedfile, and toggle whether to keep consequences (which can be really big)
filterVCF$filterVcf <- function(vcf_file, vcfout, tmpdir="./tmp", founders=NULL, bedfile=NULL, keepcsq=T, chr = NULL) 
{
    if(!is.null(chr) & !is.null(bedfile))
    {
        stop("both can't be nonnull")
    }

    chrArg = ""
    if(!is.null(chr))
    {
        chrArg = paste0(" --chr ", chr)
    }
    dir.create(tmpdir, showWarnings = F)
    ##	founderString = paste0("\n--indv ", paste(founders, collapse ="\n --indv "), "\n")
    founderString = ""
    if(!is.null(founders))
    {
	
        foundersFile = filterVCF$writeFoundersFile(tmpdir = tmpdir, founders = founders)
        founderString = paste("--keep", foundersFile)
    } 
    
    bedarg = ""
    if(!is.null(bedfile))
    {
        bedarg = paste0("--bed ", bedfile)
    } 

    csqArg = ""
    if(keepcsq)
    {
        csqArg = "--recode-INFO CSQ"
    }


    command = paste("vcftools ", "--gzvcf", vcf_file,  bedarg, chrArg, founderString, "--recode", csqArg, "-c | gzip -c >",vcfout)
    print(command)
    filterVCF$.runLoggedCommand(command, "vcflog.txt", dirname(vcfout))
}

filterVCF$writeFoundersFile <- function(tmpdir, founders)
{
    foundersFile = tempfile("afounderfile", tmpdir=tmpdir)
    write.table(founders, foundersFile, col.names=F,quote=F, row.names=F)
    return(foundersFile)
}

##filters all the vcfs in a given dir, in parallel if mc.cores>1
filterVCF$filterVcfsInDir = function(vcfDirIn, vcfDirOut,  tmpdir="./tmp", vcfExtensionIn, founders=NULL, bedfile=NULL, keepcsq=T, mc.cores=1, chr = NULL)
{
    
    dir.create(vcfDirOut, showWarnings = F)
    dir.create(tmpdir, showWarnings = F)
    vcfFiles       = stringutils$getFilesWithExtension(vcfDirIn, vcfExtensionIn)
    vcfOuts        = file.path(vcfDirOut, basename(vcfFiles))
    ii = 1:length(vcfFiles)
    func = function(i)
    {
        vcf_file = vcfFiles[i]
        vcfout   = vcfOuts[i]
        filterVCF$filterVcf(vcf_file = vcf_file, vcfout = vcfout, tmpdir = tmpdir, founders = founders, bed = bedfile, chr = chr)
    }
    mclapply(ii, FUN = func,mc.cores = mc.cores)
}

##calls a python tool which does some preprocessing of the vcf file to make it easy to read in
filterVCF$.reduceVCF = function(vcf_file, founders, vcfInfoOut=NULL, vcfPosOut=NULL, tmpdir="./tmp")
{
    if(is.null(vcfInfoOut))
    {
        vcfinfoOut   = tempfile("vcfinfo", tmpdir)
    }
    if(is.null(vcfPosOut))
    {
        vcfPosOut = tempfile("vcfpos",  tmpdir) 
    }
    foundersFile = filterVCF$writeFoundersFile(tmpdir, founders)
    command = paste("python ./genomerep/pyvcf_reduce.py", vcf_file, foundersFile, vcfInfoOut, vcfPosOut)
    filterVCF$.runLoggedCommand(command)
    return(vcfInfoOut = vcfInfoOut, vcfPosOut = vcfPosOut)
}

filterVCF$.runLoggedCommand = function(command, logfile, outputDir)
{
    
    logDir = file.path(outputDir, "logs")
    print(logDir)
    dir.create(logDir, showWarnings=F)
    ##TODO logging doesnt seem to work well with this... fix it at some point
    ##	command = paste(command, ">>&", file.path(logDir, logfile))
    print(command)
    system(command)
}
