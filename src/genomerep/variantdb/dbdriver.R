source("./genomerep/variantdb/buildVariantDb.R")
genomeData = buildGenomeData$buildAllData()
buildVariantDb$buildDbs(genomeData)
