source("./genomerep/variantdb2/buildVariantDb.R")
genomeData = buildGenomeData$buildAllData()
buildVariantDb$buildDbs(genomeData)
