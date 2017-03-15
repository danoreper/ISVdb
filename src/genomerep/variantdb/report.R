
library(data.table)
library(ggplot2)
source("./loadParams.R")
source("./genomerep/variantdb/dump_parser.R")

outputdir = fp(prop$output, "isvdb_pub")
strain.batchsize = 1
chr.batchsize    = 3
all.batchsize    = 12
    
##TODO use this in surrogatcalc
getBestAccumulator <- function(batchSize=all.batchsize, mem.gb = 6)
{
    if(prop$onCluster)
    {
        bsubCommand = bsub$get.default.killdevil.bsub(numProcessPerNode = 1, memoryLimit.GB = mem.gb, queue="day")
        accum = bsub$get.bsub.accumulator("./genomerep/variantdb/report.R", bsubCommand, batchSize=batchSize)
    } else {

        corez = prop$mnp$mc.cores
        batchSize = corez*100
        accum = bsub$get.mc.accumulator(mc.cores= corez)
    }
    return(accum)
}


## micro.analysis$getNoClusterAccum <- function()
## {
##     batchSize = prop$mnp$mc.cores*100
##     accum = bsub$get.mc.accumulator(mc.cores= prop$mnp$mc.cores)
##     return(accum)
## }

    
entropFunc <- function(strain, chr, df)
{
##    df$strain[df$strain == "A_J"] = "AJ"
    ##df$strain = gsub(df$strain, pattern = "_.*", replacement = "")

    notrans = copy(df)
    notrans$transcript_id = NULL
    notrans$transcript_name = NULL
    notrans$gene_id = NULL
    notrans$gene_name = NULL
    notrans$consequence_1 = NULL
    notrans$consequence_2 = NULL
    notrans = unique(notrans)
    ent = entropy::entropy

    entrops = notrans[,list(
        chr = chr,
        strain = strain,
        pos = pos[1],
        entrop = ent(prob)),
##        maxprob = max(prob)),
        by = "variant_id"]
    ##                            lnth = .N,
    ##      prb = paste(prob, collapse=","))
    ##, by=c("variant_id", "strain")]

    return(entrops)
}


buildColorMap <- function(avec, ncolz) 
{

    avecu = levels(avec)
    nlev = length(avecu)
    color.vec = c()

    if(ncolz%%2==1)
    {
        color.vec = as.integer(avecu)%%2
    } else {
        
        row = 0
        for(i in 1:nlev)
        {
            thecolor = (row+i)%%2
            ##flip the color if starting a new row, and each column has an even number
            if ((i-1)%%ncolz == 0)
            {
                thecolor = !thecolor
                row = row + 1
            }
            color.vec[i] = thecolor
        }
    }
    names(color.vec) = avecu
    out = color.vec[avec]
    return(out)
}


plotEntropyPerStrain <- function()
{
    print("plotting")
    outputdir = fp(prop$output,"isvdb_pub")
    dir.create(outputdir, recursive = T, showWarnings = F)

    runforstrain = function(strain)
    {
        ncolz = 4
        out = rbindlist(dump_parser$iterate(entropFunc, strain=strain, accum = getBestAccumulator(11)))
        out$chr = factor(out$chr, levels = c(1:19, "X", "Y", "MT"))
        out$color = buildColorMap(out$chr, ncolz)
        
        aplot = ggplot(out[entrop>0], aes(x=pos, y = entrop, color = color))
        ##        aplot = ggplot(out[maxprob<.99], aes(x=pos, y = maxprob, color = color))
        aplot = aplot + geom_point(size = .4)
        
        aplot = aplot + facet_wrap(~chr, ncol = ncolz)
##        aplot = aplot + ggtitle(paste0("Non-zero entropies on exons (+/-100 bp) in all chromosomes, in strain ",strain ))
        ##aplot = aplot + labs(caption = (paste0("Non-zero entropies on exons (+/-100 bp) in all chromosomes, in strain ",strain )))
        
        aplot = aplot + labs(x = "chromosomal position (bp)", y = "entropy")
        ##aplot = aplot + labs(x = "chromosomal position (bp)", y = "probability of max likelihood genotype")
        aplot = aplot + theme_bw()
        aplot = aplot + theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank())
        aplot = aplot + theme(panel.grid.major.x = element_blank(),  panel.grid.minor.x = element_blank())
        aplot = aplot + theme(panel.background = element_rect(fill = "white"))
        aplot = aplot + theme(axis.line = element_line(colour = "black"))
        aplot = aplot + theme(axis.text.x = element_text(angle = 90, hjust = 1))
        aplot = aplot + theme(strip.text = element_blank())
        aplot = aplot + theme(strip.background = element_blank())
        aplot = aplot + theme(strip.placement = "inside")
        aplot = aplot + theme(panel.spacing.y = unit(.3, "lines"))
        simpledata = unique(out[,c("chr","color"),,with=F])
        middle     = round(max(out$pos)/2)
        
        aplot = aplot + guides(fill = F, colour = F)
        aplot = aplot + geom_text(size = 6, x=middle,
                                  y = .9*max(out$entrop),
                                  ##y = .9*max(out$maxprob),
                                  hjust = .5, aes(label = chr, color = color),data = simpledata)
        ggsave(plot = aplot,
               file = fp(outputdir, paste0("entstrain_",strain,".png")),
               width = 7.5, height = 9.5 )
    }

    accum.outer = getBestAccumulator(1)
    force(runforstrain)
    accum.outer$init(func = runforstrain)
    for(strain in dump_parser$getStrains())
    {
        accum.outer$addCall(funcArgs = list(strain=strain))
    }
    accum.outer$runAll()
}

brkdownStrains <- function(strains)
{
    strainsOld = strains
    founders = c("AJ", "B6", "129", "NOD", "NZO", "CAST", "PWK", "WSB")
    
    names(founders) = c("A_J", "C57BL6J", "129S1_SvImJ", "NOD_ShiLtJ", "NZO_HlLtJ", "CAST_EiJ", "PWK_PhJ", "WSB_EiJ")    
    cclines = paste0("CC", sprintf(1:1000, fmt = "%003d"))
    names(cclines) = cclines
    all.lines = c(founders, cclines)


    strains = all.lines[as.character(strains)]
##    strains = factor(strains, levels = all.lines)

    strains = as.character(strains)

    ustrains = unique(strains)
    otherlines = ustrains[ustrains %in% founders]
    cclines    = ustrains[ustrains %in% cclines]

    strains = factor(strains, levels = all.lines)
    
    out = list(cclines = cclines, otherlines = otherlines, strains = strains)
    return(out)
}

plotEntropyPerChr <- function()
{
    print("plotting")
    outputdir = fp(prop$output,"isvdb_pub")
    dir.create(outputdir, recursive = T, showWarnings = F)
    
    runforchr = function(chr)
    {
        ncolz = 8
        out = rbindlist(dump_parser$iterate(entropFunc, chr = chr, accum = getBestAccumulator(10)))
        out$strain = factor(out$strain)
        
        strainz = brkdownStrains(out$strain)
        browser()
        out$strain = strainz$strains
                
        out$color = buildColorMap(out$strain, ncolz)

        
        aplot = ggplot(out[entrop>0], aes(x=pos, y = entrop, color = color))
        aplot = aplot + geom_point(size = .4)
        aplot = aplot + facet_wrap(~strain, ncol = ncolz)
##        aplot = aplot + ggtitle(paste0("Non-zero entropies for all strains, in exons (+/-100bp) in chromosome ",chr))
##        aplot = aplot + labs(caption = paste0("Non-zero entropies for all strains, in exons (+/-100bp) in chromosome ",chr))

        
    #    aplot = aplot + labs(x = "chromosomal position (bp)", y = "entropy")
        aplot = aplot + labs(x = "chromosomal position (bp)", y = "entropy")
        ##aplot = aplot + labs(x = "chromosomal position (bp)", y = "probability of max likelihood genotype")
        aplot = aplot + theme_bw()
        aplot = aplot + theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank())
        aplot = aplot + theme(panel.grid.major.x = element_blank(),  panel.grid.minor.x = element_blank())
        aplot = aplot + theme(panel.background = element_rect(fill = "white"))
        aplot = aplot + theme(axis.line = element_line(colour = "black"))
        aplot = aplot + theme(axis.text.x = element_text(angle = 90, hjust = 1))
        if(chr !="MT")
        {
            aplot = aplot + theme(strip.text = element_blank())
            aplot = aplot + theme(strip.background = element_blank())
        }
        aplot = aplot + theme(strip.placement = "inside")
        


#        if(chr !="MT")
#        {
            simpledata = unique(out[,c("strain","color"),,with=F])
            middle     = round(max(out$pos)/2)
                    
            aplot = aplot + geom_text(size = 4, x=middle,
                                  y = .9*max(out$entrop),
                                  ##y = .9*max(out$maxprob),
                                  hjust = .5, aes(label = strain, color = color),data = simpledata)
#        }
        
        aplot = aplot + guides(fill = F, colour = F)

        ggsave(plot = aplot,
               file = fp(outputdir, paste0("entchr_",chr,".png")),
               width = 7.5, height = 9.5 )
    }

    accum.outer = getBestAccumulator(1)
    force(runforchr)
    accum.outer$init(func = runforchr)
    for(chr in dump_parser$getChrs())
    {
##        runforchr(chr)
##        browser()
        print("done running for chr")
        accum.outer$addCall(funcArgs = list(chr=chr))
    }

    accum.outer$runAll()
}

cut2 <- function(x, breaks) {
  labels <- paste0("(",  breaks[-length(breaks)], ",", breaks[-1L], "]")
  return(factor(labels[findInterval(x, breaks)], levels=labels))
}


histfunc <- function(strain, chr, df)
{
    print(paste0(strain, "_", chr))
    df$bucket = cut2(df$prob, breaks = seq(from = .0, to = (1+.05), by =.05)-.025)
    out = data.frame(table(df$bucket))
    out$bucket =  paste0(seq(from = 0, to = 1, by = .05))
    out$strain = strain
    out$chr    = chr
    return(out)
}



plotProbHist <- function()
{
    outputdir = fp(prop$output, "isvdb_pub")
    dir.create(outputdir, recursive = T, showWarnings = F)

    ncolz = 8
    print("plotting")
    out = rbindlist(dump_parser$iterate(histfunc, accum = getBestAccumulator()))

    strainz = brkdownStrains(out$strain)
    out$strain = strainz$strains
    
    
    out$color = buildColorMap(out$strain, ncolz)
    
    print("done with iteration!!!!!!!!!!!!!!")
    
    out = out[,list(Freq = sum(Freq), color = color[1]), by = c("strain", "bucket")]
    ##    out = out[Freq>10]

   
    aplot = ggplot(out, aes(x = bucket, y=Freq,fill = color, color = color))
    aplot = aplot + geom_bar(stat="identity")
    aplot = aplot + facet_wrap(~strain, ncol = ncolz)
    ##    aplot = aplot + ggtitle("genotype probability histogram per strain")
##    aplot = aplot + labs(caption = paste0("genotype probability histogram per strain"))
    aplot = aplot + scale_y_continuous(expand = c(0,0),
           #                            limits = c(10, NA),
                                       trans = "log10")
    aplot = aplot + scale_x_discrete(name = "genotype probability", breaks=c(0, .25, .5, .75, 1))# seq(0, 1, by =.2)))

##                                       limits=c(10, max(out$Freq)))
    
    aplot = aplot + ylab("counts per probability")
##    aplot = aplot + theme_bw()
    aplot = aplot + theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank())
    aplot = aplot + theme(panel.grid.major.x = element_blank(),  panel.grid.minor.x = element_blank())
    aplot = aplot + theme(panel.background = element_rect(fill = "white"))
    aplot = aplot + theme(axis.line = element_line(colour = "black"))
    aplot = aplot + theme(axis.text.x = element_text(angle = 90, hjust = 1))
#    aplot = aplot + theme(strip.text.x = element_text(size=8))
    aplot = aplot + theme(strip.text = element_blank())


    aplot = aplot + theme(strip.background = element_blank()) ##element_rect(colour = "white", fill = "white"))
    aplot = aplot + theme(strip.placement = "inside")
##    aplot = aplot + theme(strip.text.x = element_text(aes(color = color)))
##    aplot = aplot + theme(panel.spacing.y = unit(.3, "lines"))
    middle = round(length(unique(out$bucket))/2)



    simpledata = unique(out[,c("strain","color"),,with=F])
    
    aplot = aplot + geom_text(size = 2.7, x=middle, y = log10(1.1e+06), hjust = .5, aes(label = strain, color = color),
                              data = simpledata)
    aplot = aplot + guides(fill = F, colour = F)
##    aplot = aplot + theme(axis.line = element_blank())
##    aplot = aplot + theme(panel.grid.major.y = element_blank())
##    aplot = aplot + theme(panel.grid.minor.y = element_blank())
    
    ggsave(fp(outputdir, paste0("strain_probability_hist.png")), width = 7.5, height = 9.5)
    
    print("all done")
}


tabulateConsequences <- function()
{
    tabulateConsequences.helper <- function(strain, chr, df)
    {
        
        out = data.table(data.frame(table(unlist(str_split(c(as.character(df$consequence_1),
                                                             as.character(df$consequence_2)),
                                                           pattern = "&")))))
        
        setnames(out, old = "Var1", new = "consequence")
        out$strain = strain
        out$chr = chr
        out$numrows = nrow(df)*2 ## *2 because 2 diff consequnces per row
        return(out)
    }

    out = rbindlist(dump_parser$iterate(tabulateConsequences.helper, accum = getBestAccumulator()))

    denom = unique(out[,c("strain", "chr", "numrows"),with=F])
    denom = denom[,list(numrows = sum(numrows)), by = "strain"]
    out$numrows = NULL
    setkey(denom, "strain")
    setkey(out, "strain")
    out = denom[out]
    
    
    out = out[,list(Freq = sum(Freq), numrows = numrows[1]), by = c("strain", "consequence")]
    
    
    ##    out[,fraction := Freq/sum(Freq), by = c("strain")]
    out[,fraction := Freq/(numrows), by = c("strain")]
    lows = out$fraction<.001
    highs = out$fraction>=.001
    
    lowvals  = sprintf(fmt="%.0e",  out$fraction[lows])
    highvals = sprintf(fmt="%2.3f", out$fraction[highs])
    out$fraction = as.character(out$fraction)
    out$fraction[lows] = lowvals
    out$fraction[highs] = highvals
    out$freq_fraction = paste0(out$Freq, "(", out$fraction, ")")
    fwrite(file =fp(outputdir, "full.csq.table.csv"), out)
    
##    browser()
    
    tabl = dcast(out, strain ~ consequence, value.var = "Freq")


    nonstrain.colz = which(colnames(tabl)!="strain")
    cnames = c("strain", rev(names(sort(colSums(tabl[,nonstrain.colz, with = F], na.rm = T)))))


    tabl = dcast(out, strain ~ consequence, value.var = "freq_fraction")
    tabl = tabl[,cnames, with = F]
    tabl = data.frame(tabl, check.names = F)
    tabl[is.na(tabl)]=0
    tabl = data.table(tabl)
    tabl = data.frame(tabl, check.names=F)
    
    outputdir = fp(prop$output, "isvdb_pub")
    dir.create(outputdir, showWarnings = F, recursive = T)
    write.table_with_header(
        tabl,
        sep =",",
        file= fp(outputdir, "Table_S2.csv"),
        header = "Table S2 occurences of predicted functional consequences, split by strain",
        row.names=F)
}

write.table_with_header <- function(x, file, header, ...)
{
    cat(header, '\n',  file = file)
    write.table(x, file, append = T, ...)
}
    


tabulate.het.helper <- function(strain, chr, df)
{

    notrans   = copy(df)
    notrans$consequence_1 = NULL
    notrans$consequence_2 = NULL
    notrans = unique(notrans)
    notrans$allele_1 = as.character(notrans$allele_1)
    notrans$allele_2 = as.character(notrans$allele_2)
    
    entrops = notrans[,list(maybehet = sum((allele_1!=allele_2)*prob)>=.248),
                      by=c("variant_id")]
    
    out = data.frame(numhet = sum(entrops$maybehet), numhom = sum(!entrops$maybehet))
    out$strain = strain
    out$chr = chr
    return(out)
}


tabulateHeterozygosity <- function(pct.thresh)
{

    out = rbindlist(dump_parser$iterate(tabulate.het.helper, accum = getBestAccumulator()))
    browser()
    ##hetzyg = fread("../output/isvdb_pub/strain.heterozygosityTable.csv")
    
    out$fraction_hom= out$numhet/(out$numhet+out$numhom)
    out$fraction_hom= as.numeric(sprintf(out$fraction, fmt = "%2.3f"))
    out = out[,c(c("strain","chr"), setdiff(colnames(out), c("strain","chr"))),with = F]
    
    setnames(out, old="fraction_hom", new = "pct_het")
    outputdir = fp(prop$output, "isvdb_pub")
    dir.create(outputdir, showWarnings = F, recursive = T)
    data.table::fwrite(out, fp(outputdir, "strain.heterozygosityTable.csv"))

    strainz = brkdownStrains(out$strain)
    out = out[strain %in% strainz$cclines]
    
    totals = out[ ,list(numhet = sum(numhet), numhom = sum(numhom)), by = strain]
    totals$all_chr = totals$numhet/(totals$numhet+totals$numhom)
    totals$numhet = NULL
    totals$numhom = NULL
    setkey(totals, "strain")

    smry = dcast(out, strain~chr, value.var = "pct_het")
    setkey(smry, "strain")
    bigtable = smry[totals]
    rownames(bigtable) = bigtable$strain
    bigtable$strain = NULL
    bigtable$all_chr = sprintf("%.3f", bigtable$all_chr)
    
    fle = fp(outputdir, "Table_S1.csv")
    ##    xtable::print.xtable(xtable::xtable(bigtable, caption = 'Table S2 Heterozygosity fraction per strain', caption.placement = 'top'), file = fle)
    write.table_with_header(bigtable, fle, header="Table S1 residual heterozygosity fraction per strain", sep=",", row.names = F)
    
    ## setnames(totals, old = "all_chr", new = "pct_het")
    ## fle = fp(outputdir, "Table_S1.tex")
    ## totals$pct_het= sprintf("%.3f", totals$pct_het)
    ## xtable::print.xtable(xtable::xtable(totals, caption = 'Table S1 Heterozygosity fraction per strain', caption.placement = 'top'), file = fle)
}


tabulate.het.helper <- function(strain, chr, df)
{
    write.table_with_header <- function(x, file, header, ...)
    {
        cat(header, '\n',  file = file)
        fwrite(x, file, append = T, ...)
    }

    
    
    notrans               = copy(df)
    notrans$consequence_1 = NULL
    notrans$consequence_2 = NULL
    notrans               = unique(notrans)
    notrans$allele_1      = as.character(notrans$allele_1)
    notrans$allele_2      = as.character(notrans$allele_2)
    
    entrops = notrans[,list(maybehet = sum((allele_1!=allele_2)*prob)>=.248),
                      by=c("variant_id")]
    
    out = data.frame(numhet = sum(entrops$maybehet), numhom = sum(!entrops$maybehet))
    out$strain = strain
    out$chr = chr
    return(out)
}



fixfiles <- function(
                     strains = dump_parser$getStrains(),
                     chrs    = dump_parser$getChrs(strain = strains[1]))


{

    fixfunc = function(strain, chr, df)
    {
        col.names = c("variant_id", "pos", "allele_1", "allele_2", "prob", "is_max", "consequence_1", "consequence_2", "gene_name", "transcript_name")
        headerlines = paste0("strain:", strain, ",chr:",chr, "\n")
        headerlines = paste0(headerlines, paste0(paste(col.names, collapse=","),"\n"))
        dump_parser$write(strain, chr, df, headerlines)
        
        
    }
    dump_parser$iterate(fixfunc, accum = getBestAccumulator(20))
}





fixdiplos <- function(
                     strains = dump_parser$getStrains(),
                     chrs    = dump_parser$getChrs(strain = strains[1]))


{

    fixfunc = function(strain, chr, df)
    {
        col.names = c("variant_id", "pos", "founder_1", "founder_2", "prob", "gene_name")
        headerlines = paste0("strain:", strain, ",chr:",chr, "\n")
        headerlines = paste0(headerlines, paste0(paste(col.names, collapse=","),"\n"))
        dump_parser$write(strain, chr, df, headerlines)
    }
    dump_parser$iterate(fixfunc, accum = getBestAccumulator(20))
}

