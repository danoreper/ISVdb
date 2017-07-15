## This code can be modified to get diplotypes, inbred lines themselves rather than crosses, and can beeasily run in parallel if you need it to (e.g.., to rapidly generate lots of random crosses)
## Let me know if you want to do those things.


source("./genomerep/variantdb2/dump_parser.R")
root = "/nas/depts/006/valdar-lab/PUBLIC/isvdb/"

##dbdir = fp(root, "exon_1410") ##MGP 1410, Just the exons-- manipulating the exon data is much faster than whole genome.
##dbdir = fp(root, "full_1504") ##MGP 1504, Whole genome
dbdir = fp(root, "exon_1504")

db    = db_builder$get.db.lite(dbdir)

allchrs = db$genotype$getChrs()


strain1s = c("CC001", "CC002", "CC003")
strain2s = c("CC004", "CC005", "CC006")
chr = "7"

crosses = list()
for(i in 1:length(strain1s))
{
    print(paste0("cross", i))
    strain1 = strain1s[i]
    strain2 = strain2s[i]

    cross = db$read(strain1 = strain1, strain2 = strain2, type = "diplotype", phased = T, chr = chr)
    cross_wide = db_builder$toWideMappingFormat(cross)
    crosses[[i]] = cross_wide
}
    
    


