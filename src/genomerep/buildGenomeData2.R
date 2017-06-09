# TODO: Add comment
# 
# Author: doreper
###############################################################################
library("GenomicRanges")
library("IRanges")
library("data.table")
library("Biostrings")
library("rtracklayer")

source("./genomerep/buildImprintedGenes.R")
source("./loadParams.R")
source("./utils.R")

buildGenomeData = new.env()

buildGenomeData$buildAllData=function(buildVariantDb = F)
{
    ##	source("./strainselection/hetcalc/buildIBD.R")
    ##	buildGenomeData$ibdData  = buildIBD$build(newickDir = newickDir)

    dir.create(prop$output)
    dir.create(prop$tmpdir)
    genomeData = new.env(hash=T)

    genomeData$founders = read.table(dat(prop$genome$foundersMap),header=T, sep=",")$founder
    ##genomeData$karyotype = buildGenomeData$getKaryotype(dat( prop$genome$dnaReferenceFile))
    genomeData$karyotype = buildGenomeData$getKaryotype(dat( prop$genome$karyotype))

    genomeData$snords = buildGenomeData$getSnords()
    
    if(prop$genome$rebuildImprinted) 
    {
        ##writen to imprinted file location
        buildGenomeData$rebuildImprintedFile()
    }
    genomeData$imprintedGenes = buildGenomeData$loadImprintedGranges(dat(prop$genome$imprintedGenesFile))

    ##TODO wrap up the slow part
    print("getting exons")
    if(prop$genome$rebuildImprintedTrancriptsFile | prop$genome$loadexons)
    {
        print(dat(prop$genome$exonReferenceFile))
##        genomeData$exonFeatures               = import(dat( prop$genome$exonReferenceFile), "gtf")
        genomeData$exons                      = buildGenomeData$getExons(dat(prop$genome$exonReferenceFile))

        print("done getting exons")

        imprintedGenes = as.data.frame(unique(genomeData$imprintedGenes))
        genomeData$imprintedExonFeatures = genomeData$exons[genomeData$exons$gene_id %in% imprintedGenes$ensembl_gene_id]

        ##transcript regions rather than actual transcript names... find out from Fernando if this is what we want
        ## genomeData$imprintedExonFeatures      = genomeData$exons[subjectHits(findOverlaps(genomeData$imprintedGenes, genomeData$exons))]
        
        genomeData$imprintedGeneTranscriptIDs = unique(genomeData$imprintedExonFeatures$transcript_id)
        write.table(genomeData$imprintedGeneTranscriptIDs, file=dat(prop$genome$imprintedTranscriptsFile), row.names=F, quote=F, col.names=F)
        genomeData$exons=NULL
    }
    print("done with slow")
    
    if(buildVariantDb)
    {
        buildGenomeData$buildVariantDbs(genomeData)
    }
    return(genomeData)
}

buildGenomeData$getSnords <- function()
{
    snords = dat(fp(prop$genome$snord))
    dir(snords, "*.csv")
    dfs = list()
    for(snord in dir(snords,"*.csv"))
    {
        id  = sub("\\.csv","",snord )
        df = fread(fp(snords,snord))
        df$snord = id
        dfs = util$appendToList(dfs, df)
    }
    df = rbindlist(dfs)
    return(df)
}

buildGenomeData$rebuildImprintedFile <- function()
{
    outfile         = dat( prop$genome$imprintedGenesFile)
    fullImprintInfo = buildImprintedGenes$collateImprintInfo()
    write.table(fullImprintInfo, outfile, sep="\t", row.names = F)
}

#Returns an IRanges object containing metadata as well as all genes in imprinted genes file
buildGenomeData$loadImprintedGranges <- function(imprintedGenesFile) 
{
    geneMetaData = read.table(imprintedGenesFile, header=T, sep="\t")
    strand = as.character(geneMetaData$strand)
    strand[geneMetaData$strand=="1"]   = "+"
    strand[geneMetaData$strand=="-1"]  = "-"
    
    gr = GRanges(seqnames = geneMetaData$chromosome_name, 
                 ranges = IRanges(geneMetaData$start_position, geneMetaData$end_position, names = geneMetaData$mgi_symbol),
                 strand=strand,
                  geneMetaData[,!(colnames(geneMetaData) %in% c("chromosome_name", 
                                                                       "start_position",
                                                                       "end_position",
                                                                       "strand"))]
                 )
    return(gr)
}

buildGenomeData$buildGr <- function(start, end, strand, chr, names = chr)
{
    ir = IRanges(start = start,
                 end   = end,
                 names = names)

    gr = GRanges(ranges=ir,
                 strand   = strand,
                 seqnames = chr)

    return(gr)
}
    

buildGenomeData$getExons <- function(gtfFile = dat( prop$genome$exonReferenceFile)) 
{
    agtf=rtracklayer::import(gtfFile)
    print("imported exon")
    agtf$score = NULL #this is NA anyway
    exons = agtf[agtf$type=="exon",]
    return(exons)
}

buildGenomeData$getGeneInfo <- function(gtfFile = dat( prop$genome$exonReferenceFile))
{
    print("getting exons")
    exonz = buildGenomeData$getExons(gtfFile)
    exonz = data.table(as.data.frame(exonz))

    print("breaking down by gene")
    geneRanges = exonz[ ,list(start=min(start),
                              end = max(end),
                              chrom = paste(unique(seqnames), collapse=","),
                              strand = paste(unique(strand), collapse=","),
                              gene_name = paste(unique(gene_name), collapse=",")),
                       by = gene_id]

    print("keying by gene_id")
    setkey(geneRanges, "gene_id")
    return(geneRanges)
}

## ##create a data frame containing length per chromosome
## buildGenomeData$getKaryotype <- function(referenceFile =  dat(prop$genome$dnaReferenceFile)) 
## {
##     refstring = readDNAStringSet(referenceFile)
##     chrname   = strsplit(names(refstring), " ")
##     chrname   = unlist(lapply(chrname, "[", 1))
##     lengths   = width(refstring)
##     rm(refstring)
##     gc()
##     karyotypeFrame = data.frame(chrname = chrname , len=lengths)
##     rownames(karyotypeFrame) = chrname
##     return(karyotypeFrame)
## }
buildGenomeData$getKaryotype <- function(referenceFile =  dat(prop$genome$karyotype)) 
{
    ##    browser()
    print(referenceFile)
    a = fread(referenceFile)
    setnames(a, old = c("V1", "V2"), new= c("chrname", "len"))
    a = a[chrname %in% c("chrX", "chrY", "chrM", paste0("chr",1:19))]
    a$chrname = gsub(a$chrname, pattern = "chr", replacement = "")
    a$chrname[a$chrname =="M"] = "MT"
    setkey(a, "chrname")
    rownames(a) = a$chrname
##    browser()
    return(a)
}

buildGenomeData$writeToBed <- function(ranges, exonBedFile)
{
    export(ranges, exonBedFile,format="bed")	
}


