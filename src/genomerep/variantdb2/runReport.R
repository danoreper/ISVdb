source("./genomerep/variantdb2/report.R")
print("ON CLUSTER????")
print(prop$onCluster)




##TODO uncomment this
db = db_builder$get.db.lite(dat(prop$variantdb$dbfolder))
##db = db_builder$get.db.lite(dat(fp("isvdb", "exon_1410")))

pracma::tic()
db.report$tabulateConsequences(db)
pracma::toc()



pracma::tic()
db.report$tabulateHeterozygosity(db, .25)
pracma::toc()


pracma::tic()
counter = db.report$plotEntropyPerChr(db)
pracma::toc()


db.report$buildLatex(db)


pracma::tic()
counter = db.report$plotEntropyPerStrain(db)
pracma::toc()



pracma::tic()
counter = db.report$plotProbHist(db)
pracma::toc()


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



