## Given a cc probs file of HMM output specifying diplotypes at every MegaMUGA marker,
## lets the user query a set of positions between those markers for a new lineaerly interpolated
## set of diplotype probabilities
##
## To use, build a probData object, and then query it repeatedly, or all at once with a vector of positions.
##
## Author: doreper
###############################################################################
library("igraph")
library("IRanges")
library("stringr")
library("data.table")
library("GenomicRanges")

source("./stringutils.R")

buildFounderProbs = new.env(hash=T)

##probData is the data structure that is the result of buildFounderProbs$build
##variant_id: a vector of labels for the queried positions.
##chrom: a vector of chromosomes
##pos : an equal length vector of positions
##TODO: include a start and end rather than just a position for vairant that span some range
buildFounderProbs$hapProbPerPosition <- function(variant_id, chrom, pos, probData) 
{
    if(length(unique(variant_id))!=length(pos))
    {
        
        stop("bad variant Ids")
    }
    ## build granges for the queried positions
    variantIR           = GRanges(ranges = IRanges(pos, pos), seqnames=chrom, variant_id = variant_id) 

    ##overlap the queried positions with the left-right pairs of markers
    overlapsPerHap      = findOverlaps(variantIR, probData$ir)
    hapIndex            = subjectHits(overlapsPerHap)   # an index into the haplotype range that overlaps with each variant
    varIndex            = queryHits(overlapsPerHap)
    varPositions        = start(variantIR[varIndex])
    ##the markers from the haplotype range information that corresponds to each overlap. 
    ##Markers come with the name of the marker, which can be used to loop up the appropriate marker probs
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

##Get the probs file for line out of the dipolotype probs file dir.
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

## build a probData object for querying specific positions later. Build this once.
## and then call hapProbPerPosition using it repeatedly.
##probsfile: a diplotype file for a single strain
## founders: a vector of the full founder strain names (A-H).
## karyotype: a datatable in which the key is the chromosome name, and the value is a length per chromosome
buildFounderProbs$build <- function(probsFile, founders, karyotype, chr = NULL) 
{
    markerProbs = buildFounderProbs$.buildMarkerProbs(probsFile = probsFile, founders) 
    founderPairs = markerProbs$founderPairs
    markerProbs  = markerProbs$markerProbs
    print("building prob ranges")
    ir          = buildFounderProbs$.buildProbRanges(markerProbs, karyotype, chr)
    setkey(markerProbs, "Marker.Name")

    print("returning")
    return(list(markerProbs=markerProbs, ir=ir, founderPairs = founderPairs))
}

##build a data frame indexed by marker name with a set of diplotype probabilities per marker name
buildFounderProbs$.buildMarkerProbs <- function(probsFile, founders) 
{
    ## strainNameAndDescent = buildFounderProbs$parseProbFilename(probsFile)
    ## print(strainNameAndDescent)
    probFrame = try((fread(probsFile, header=T, sep=",", stringsAsFactors = F)))

    old.and.new = buildFounderProbs$.renameDiploCols(probFrame,founders)
    setnames(probFrame, old = c("marker",      "position(B38)", "chromosome", old.and.new[["oldnames"]]),
                        new = c("Marker.Name", "Pos",           "Chrom",      old.and.new[["newnames"]]))
    founderPairs = old.and.new[["newnames"]]
    probFrame = na.omit(probFrame)
    probFrame$Chrom[probFrame$Chrom=="M"] = "MT"

    probFrame = data.table(probFrame, key=c("Chrom"))
    gc()
    return(list(markerProbs = probFrame, founderPairs = founderPairs))
}

##build up a data frame with full founder strain names rather than just letters, for ease of debugging.
buildFounderProbs$.renameDiploCols <- function(probFrame, founders)
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


##builds ranges with known left marker and right marker probabilities, which will be used to identify in between which two markers a queried position lands.
buildFounderProbs$.buildProbRanges <- function(markerProbs, karyotype, chrs = unique(markerProbs$Chrom)) 
{
    if(is.null(chrs))
    {
        chrs =  unique(markerProbs$Chrom)
    }
    grs = list()
    counter     = 1
    
    for(chr in chrs)
    {
        grs[[counter]] = buildFounderProbs$.buildProbRange(markerProbs = markerProbs, chr = chr, karyotype = karyotype)
        counter = counter + 1
    }
    grs = do.call(c, grs)
    return(grs)
}

##Given a set of markerProbs on a single chromosome, builds a grange for every pair of markers. Looks up chromosome length in karyotype
buildFounderProbs$.buildProbRange <- function(markerProbs, chr, karyotype) 
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
    
    karyoLen = karyotype[chr]$len
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


##Get the strain names and their descent out of a probs file
buildFounderProbs$parseProbFilename <- function(probFileName)
{
    strainNameAndDescent = unlist(strsplit(basename(probFileName),split=".csv"))
    strainNameAndDescent = (strsplit(strainNameAndDescent, split = "_"))
    strainName = lapply(strainNameAndDescent, "[", 1)
    descent = lapply(strainNameAndDescent, "[", 2)
    ret = list(strainName = strainName, descent = descent)
    return(ret)
}
