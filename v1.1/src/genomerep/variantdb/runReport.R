source("./genomerep/variantdb/report.R")

print("ON CLUSTER????")
print(prop$onCluster)


## ##outputdir = fp(prop$output, "isvdb_pub")
buildLatex <- function()
{
    sink(fp(outputdir, "progchr.tex"))
    ##chrs = c(1:19, "X", "Y", "MT")
    
    chrs = dump_parser$getChrs()
    for(chr in chrs)
    {
        astring = 
        paste0("
\\begin{figure*}[htbp]
\\renewcommand{\\familydefault}{\\sfdefault}\\normalfont
\\centering
\\includegraphics[width= \\textwidth]{entchr_",chr,"}
\\caption[chr ", chr,", non-zero entropies] {chr ", chr,", non-zero entropies in exons (+/-100bp) in all strains. Each point
corresponds to the entropy of a variant at that position along the chromosome}
\\end{figure*}
\\clearpage
")      
        cat(astring)
    }
    sink()

    
    
    sink(fp(outputdir, "progstrain.tex"))
    strains = dump_parser$getStrains()
    strain_safes = brkdownStrains(strains)
    for(i in 1:length(strains))
    {
        strain = strains[i]
        ##        strain_safe = gsub(strain, pattern = "_", replacement = "\\\\_")
        strain_safe = strain_safes$strains[i]
        astring = 
        paste0("
\\begin{figure*}[htbp]
\\renewcommand{\\familydefault}{\\sfdefault}\\normalfont
\\centering
\\includegraphics[width= \\textwidth]{entstrain_",strain,"}
\\caption[strain ", strain_safe,", non-zero entropies] {strain ", strain_safe,", non-zero entropies in exons (+/-100 bp) in all chromosomes. Each point corresponds to the entropy of a variant at that position along the chromosome } 
\\end{figure*}
\\clearpage
")
        cat(astring)
    }
    sink()
}

pracma::tic()
tabulateConsequences()
pracma::toc()
browser()


pracma::tic()
tabulateHeterozygosity(.25)
pracma::toc()

browser()

pracma::tic()
counter = plotEntropyPerChr()
pracma::toc()
browser()


buildLatex()
browser()


browser()
pracma::tic()
counter = plotEntropyPerStrain()
pracma::toc()

browser()



pracma::tic()
counter = plotProbHist()
pracma::toc()


## \begin{minipage}{\linewidth}
##     \centering%
##     \includegraphics[width=1\textwidth]{{/Users/wvaldar/Dropbox/PsyDiallel/Weight/pre2/TwoPlot <- weight <- pre.png}}
##     \figcaption{\label{}Body Weight observed and predicted phenotype values }
##     \end{minipage}

##    \includegraphics[width=1\textwidth]{{/Users/wvaldar/Dropbox/PsyDiallel/Weight/pre2/TwoPlot <- weight <- pre.png}}

##outputdir = fp(prop$output, "isvdb_pub")

## buildLatex <- function()
## {
##     sink(fp(outputdir, "progchr.tex"))
##     ##chrs = c(1:19, "X", "Y", "MT")
    
##     chrs = dump_parser$getChrs()
##     for(chr in chrs)
##     {
##         astring = 
##             paste0("
## \\begin{minipage}{\\linewidth}
## \\centering%
## \\includegraphics[width= \\textwidth]{entchr_",chr,"}
## \\figcaption[chr ", chr,", non-zero entropies] {chr ", chr,", non-zero entropies in exons (+/-100bp) in all strains. Each point
## corresponds to the entropy of a variant at that position along the chromosome}
## \\end{minipage}
## \\clearpage
## ")      
##         cat(astring)
##     }
##     sink()
    
    
##     sink(fp(outputdir, "progstrain.tex"))

##     strains = dump_parser$getStrains()

##     for(strain in strains)
##     {
##         strain_safe = gsub(strain, pattern = "_", replacement = "\\\\_")
##         astring = 
##         paste0("
## \\begin{minipage}{\\linewidth}
## \\centering%
## \\includegraphics[width= \\textwidth]{entstrain_",strain,"}
## \\figcaption[strain ", strain_safe,", non-zero entropies] {strain ", strain_safe,", non-zero entropies in exons (+/-100 bp) in all chromosomes. Each point corresponds to the entropy of a variant at that position along the chromosome } 
## \\end{minipage}
## \\clearpage
## ")
##         cat(astring)
##     }
##     sink()
## }

## buildLatex()


## > head(strains)


hetzyg = fread(fp(outputdir, "strain.heterozygosityTable.csv"))
cc.zyg = hetzyg[grepl(hetzyg$strain, pattern="CC")]

pcthet.cc = sum(cc.zyg$numhet)/(sum(cc.zyg$numhet)+sum(cc.zyg$numhom))
print(paste0("pct.cchet: ",pcthet.cc))




csq = fread(fp(outputdir, "full.csq.table.csv"))
csq.cc = csq[grepl(csq$strain, pattern = "CC")]
denom = sum(csq.cc$Freq)
print(paste0("pct cc reference:", print(sum(csq.cc$Freq[csq.cc$consequence == "reference"])/denom)))
print(paste0("pct cc intron_variant:", print(sum(csq.cc$Freq[csq.cc$consequence == "intron_variant"])/denom)))
print(paste0("pct cc downstream_gene_variant:", print(sum(csq.cc$Freq[csq.cc$consequence == "downstream_gene_variant"])/denom)))
print(paste0("pct cc nc_transcript_variant:", print(sum(csq.cc$Freq[csq.cc$consequence == "nc_transcript_variant"])/denom)))
print(paste0("pct cc upstream_gene_variant:", print(sum(csq.cc$Freq[csq.cc$consequence == "upstream_gene_variant"])/denom)))
print(paste0("pct cc stop_gained:", print(sum(csq.cc$Freq[csq.cc$consequence == "stop_gained"])/denom)))
print(paste0("pct cc stop_lost:", print(sum(csq.cc$Freq[csq.cc$consequence == "stop_lost"])/denom)))

