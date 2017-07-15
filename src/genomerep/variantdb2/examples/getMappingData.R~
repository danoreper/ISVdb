## This code can be modified to get diplotypes, inbred lines themselves rather than crosses, and can beeasily run in parallel if you need it to (e.g.., to rapidly generate lots of random crosses)
## Let me know if you want to do those things.


source("./genomerep/variantdb2/dump_parser.R")
root = "/nas/depts/006/valdar-lab/PUBLIC/isvdb/"

##dbdir = fp(root, "exon_1410") ##MGP 1410, Just the exons-- manipulating the exon data is much faster than whole genome.
##dbdir = fp(root, "full_1504") ##MGP 1504, Whole genome
dbdir = fp(root, "exon_1504")

db    = db_builder$get.db.lite(dbdir)
cross = list()
for(chr in db$genotype$getChrs())
{
    cross[[chr]] = db$read(strain1 = "NOD_ShiLtJ", strain2 = "CC001", type = "genotype", phased = T, chr = chr)
}

cross = do.call(rbind, cross)
