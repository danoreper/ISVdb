                                        # Given a cc probs file of the form that ping generates, that is a readout of HMM probs per megamuga marker, build a structure that 
                                        # can be efficiently queried for the diplotype probability at given chromosome and position
                                        # 
                                        # Author: doreper
###############################################################################
library("igraph")
library("IRanges")
library("stringr")
library("data.table")
library("GenomicRanges")

source("./stringutils.R")

buildFounderProbs = new.env(hash=T)

buildFounderProbs$renameDiploCols <- function(probFrame, founders)
{
    oldColNames = expand.grid(LETTERS[1:8], LETTERS[1:8])
    oldColNames = paste0(oldColNames$Var1, oldColNames$Var2)
    goodNames = (oldColNames %in% colnames(probFrame))
    oldColNames = oldColNames[goodNames]
    
    newColNames = expand.grid(founders, founders)
    newColNames = paste0(newColNames$Var1,".", newColNames$Var2)
    newColNames = newColNames[goodNames]
    
    old.and.new = list(oldnames = oldColNames, newnames = newColNames)
    return(old.and.new) 
}

buildFounderProbs$buildMarkerProbs <- function(probsFile, founders) 
{
    strainNameAndDescent = buildFounderProbs$parseProbFilename(probsFile)
    print(strainNameAndDescent)
    probFrame = try((fread(probsFile, header=T, sep=",", stringsAsFactors = F)))

    old.and.new = buildFounderProbs$renameDiploCols(probFrame,founders)
    setnames(probFrame, old = c("marker",      "position(B38)", "chromosome", old.and.new[["oldnames"]]),
                        new = c("Marker.Name", "Pos",           "Chrom",      old.and.new[["newnames"]]))
    founderPairs = old.and.new[["newnames"]]
    probFrame = na.omit(probFrame)
    probFrame$Chrom[probFrame$Chrom=="M"] = "MT"

    probFrame = data.table(probFrame, key=c("Chrom"))
    gc()
    return(list(markerProbs = probFrame, founderPairs = founderPairs))
}

##builds ranges with known left marker and right marker probabilities, which will be used to interpolate probabilities in between when queried
buildFounderProbs$buildProbRanges <- function(markerProbs, karyotype, chrs = unique(markerProbs$Chrom)) 
{
    ##
    grs = list()
    counter     = 1
    
    for(chr in chrs)
    {
        grs[[counter]] = buildFounderProbs$buildProbRange(markerProbs = markerProbs, chr = chr, karyotype = karyotype)
        counter = counter + 1
    }
    grs = do.call(c, grs)
    return(grs)
}

buildFounderProbs$buildProbRange <- function(markerProbs, chr, karyotype) 
{
    subTable = markerProbs[J(Chrom=chr)]
    if(length(table(table(subTable$Marker.Name)))>1)
    {
        stop("bad markers, at least one marker has a redundant name")
    }
    if(!(chr %in% rownames(karyotype)))
    {
        
        stop(paste0("missing karyotype for ", chr))
    } 
    
    setkey(subTable, "Pos") #Orders the data frame by position
    
    karyoLen = karyotype[chr,"len"]
    numFenceposts = nrow(subTable)+ 1
    starts      = c(0, subTable$Pos)
    ends        = c(subTable$Pos,   max(karyoLen, subTable$Pos[nrow(subTable)]))
    gr = GRanges(ranges = IRanges(starts, ends), 
                 seqnames = chr, 
                 leftMarker  = c(subTable$Marker.Name[1], as.character(subTable$Marker.Name)), 
                 rightMarker = c(as.character(subTable$Marker.Name), subTable$Marker.Name[nrow(subTable)]))
    print("built Ranges")
                                        #	
    return(gr)
}

buildFounderProbs$build <- function(probsFile, founders, karyotype) 
{
                                        #	Rprof(line.profiling=T)
    markerProbs = buildFounderProbs$buildMarkerProbs(probsFile = probsFile, founders) 
                                        #	Rprof(NULL)
    founderPairs = markerProbs$founderPairs
    markerProbs  = markerProbs$markerProbs
    print("building prob ranges")
    ir          = buildFounderProbs$buildProbRanges(markerProbs, karyotype)
    setkey(markerProbs, "Marker.Name")

    print("returning")
    return(list(markerProbs=markerProbs, ir=ir, founderPairs = founderPairs))
}


##probData is the data structure that is the result of buildFounderProbs$build 
##POS column, and CHROM column
##TODO: include a start and end rather than just a position for vairant that span some range
buildFounderProbs$hapProbPerPosition <- function(variant_id, chrom, pos, probData) 
{
    if(length(unique(variant_id))!=length(chrom))
    {
        
        stop("bad variant Ids")
    }
    variantIR           = GRanges(ranges = IRanges(pos, pos), seqnames=chrom, variant_id = variant_id) 
    
    overlapsPerHap      = findOverlaps(variantIR, probData$ir)
    hapIndex            = subjectHits(overlapsPerHap)   # an index into the haplotype range that overlaps with each variant
    varIndex            = queryHits(overlapsPerHap)
    varPositions        = start(variantIR[varIndex])
                                        #the markers from the haplotype range information that corresponds to each overlap. 
                                        #Markers come with information about which strain they are measured in, as well as the name of the marker
    boundingMarkers     = probData$ir[hapIndex]

    leftProbData    = probData$markerProbs[J(Marker.Name=boundingMarkers$leftMarker)]
    rightProbData   = probData$markerProbs[J(Marker.Name=boundingMarkers$rightMarker)]
    
    leftDistribution       = leftProbData[,probData$founderPairs, with = F]
    rightDistribution      = rightProbData[,probData$founderPairs, with = F]
    
    slope             = (rightDistribution-leftDistribution)/(rightProbData$Pos - leftProbData$Pos)
    slope[boundingMarkers$leftMarker==boundingMarkers$rightMarker] = 0; #special code to handle 0 to first marker, and last marker to karyotype length
    
    varProbDist       = leftDistribution + slope* (varPositions - leftProbData$Pos)

    varProbDist         = cbind(variant_id=variantIR[varIndex]$variant_id, varProbDist)
    varProbDist         = varProbDist[!duplicated(varProbDist$variant_id),]
    return(varProbDist)
}

buildFounderProbs$parseProbFilename <- function(probFileName)
{
    strainNameAndDescent = unlist(strsplit(basename(probFileName),split=".csv"))
    strainNameAndDescent = (strsplit(strainNameAndDescent, split = "_"))
    strainName = lapply(strainNameAndDescent, "[", 1)
    descent = lapply(strainNameAndDescent, "[", 2)
    ret = list(strainName = strainName, descent = descent)
    return(ret)
}

buildFounderProbs$getProbsFileForLine <- function(line, rilHaplotypeProbsDir= dat(prop$genome$rilHaplotypeProbsDir))
{
    line = as.character(line)
    line = strsplit(line, "/")
    line = unlist(lapply(line, "[", 1))
    probsFiles = stringutils$getFilesWithExtension(rilHaplotypeProbsDir, ".csv")
    probsFile = probsFiles[grepl(line, probsFiles)]
    if(length(probsFile)!=1)
    {
        stop(paste("matching probs files:", paste(probsFile, collapse = ",")))
    }
    return(probsFile)
}
