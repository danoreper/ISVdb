# Builds data frame of imprinted gene info, with names and positions taken from ensembl if possible. We use a litfile, including those genes
# that were specified in mousebook, and a file specifying those gene found to be brain expressed in the crowley paper
# Author: doreper
###############################################################################

library("biomaRt")
library(data.table)

buildImprintedGenes = new.env(hash=T)

#builds a data frame containing all the imprinted gene information
buildImprintedGenes$collateImprintInfo <- function(
  imprintFile     = dat( prop$genome$imprintCrowleyExperimentFile),
  litFile         = dat( prop$genome$imprintLitFile),
  snord115.file   = dat( prop$genome$snord115.file),
  snord116.file   = dat( prop$genome$snord116.file),
  ensemblhost     = prop$genome$biomartHostForDNAReference) 
{
    ##set up biomart access
    listMarts(host= ensemblhost) 
    ensembl=useMart(host = ensemblhost, "ENSEMBL_MART_ENSEMBL") #corresponds to 38.75 release, the name string might have to change for other releases
    listDatasets(ensembl)
    ensembl = useDataset("mmusculus_gene_ensembl",mart=ensembl)
    filters = listFilters(ensembl)
    attributes = listAttributes(ensembl)
    

    ##the crowley supplementary table imprinted measured genes
    imprinted =  data.table(read.table(imprintFile,header=T,sep=",", skip=1),  key="Ensembl.name")
    
    literature = data.table(read.table(litFile,header=T,sep=","),key="Ensembl.Gene.ID")
    litOnly    = setdiff(literature$Ensembl.Gene.ID, imprinted$Ensembl.name)
    litOnly    = literature[litOnly]
    nLit       = nrow(litOnly)
    nImprint   = nrow(imprinted)
    
#collate the supplementary table imprinted genes, and the literature known imprinted genes into a single data frame
#we collect the bare minimum of what is necessary from these tables, and look up the rest from ensembl mart, as that is an online resource that is more current.
    impGenes =data.frame(
      ensembl_gene_id = c(as.character(litOnly$Ensembl.Gene.ID), as.character(imprinted$Ensembl.name)), 
      mgi_symbol      = c(as.character(litOnly$Gene), as.character(imprinted$Gene)),
      brainImprinted  = c(rep(F,nLit), rep(T, nImprint)),
      strainEffect    = c(rep(F,nLit), as.character(imprinted$Strain.effect.))=="Yes",
      Expressed.allele = c(as.character(litOnly$Expressed.allele), as.character(imprinted$Expressed.allele)))

    impGenes = data.table(impGenes,key="ensembl_gene_id")
    impGenes = impGenes[!mgi_symbol %in% c("Snord115","Snord116")]
    
    
    ##enMart  = getGene(ensembl.gene.ID,type="ensembl_gene_id", mart=ensembl) 
    enMart = getBM(attributes=c("ensembl_gene_id",
                                "mgi_symbol",
                                "description",
                                "chromosome_name",
                                "strand",           
                                "start_position",
                                "end_position"),
                   filters="ensembl_gene_id", values=impGenes$ensembl_gene_id,mart=ensembl)

    enMart = data.table(enMart, key="ensembl_gene_id")
    fullInfo                 = impGenes[enMart]
    ##use the ensembl mart version of IDs after joining the two tables, as it is more current
    fullInfo$mgi_symbol      = fullInfo$i.mgi_symbol 
    fullInfo$i.mgi_symbol      = NULL
	
    
    ##some number of ensembl ids wont exist in ensembl mart, lets examine them
    allEnsembl   = unique(c(as.character(imprinted$Ensembl.name), as.character(literature$Ensembl.Gene.ID)))
    missingEnsembl = setdiff(allEnsembl, fullInfo$ensembl_gene_id)
    print("following ensembl ids were not found in biomart... they are deprecated as of my last manual check against ensembl page")
    print(missingEnsembl)


    ##This has to be done separately since the IDS are all bogus for snord115 and 116
    
    snord115 = read.table(snord115.file, sep=",", nrows=0, header=T)
    snord115$snordtype = "snord115"

    snord116 = read.table(snord116.file, sep=",", nrows=0, header=T)
    snord116$snordtype = "snord116"

    snords = rbind(snord115, snord116)
    snords$location = as.character(snords$location)
    
    pos=  unlist(lapply(strsplit(snords$location,split=":"), "[",2))
    start = unlist(lapply(strsplit(pos, split = "-"), "[", 1))
    end = unlist(lapply(strsplit(pos, split = "-"), "[", 2))

    snords = data.frame(ensembl_gene_id = snords$id_with_url,
                        mgi_symbol      = paste0(snords$name),
                        brainImprinted  = F,
                        strainEffect    = F,
                        snordType       = snords$snordtype,
                          description     = snords$snordtype,
                          chromosome_name = unlist(lapply(strsplit(snords$location,split=":"), "[",1)),
                          strand           =  unlist(lapply(strsplit(snords$location,split=":"), "[",3)),
                          start_position   = start,
                          end_position     = end,
                          Expressed.allele = NA)

    fullInfo = rbind(fullInfo, snords, fill = T)
    #redundant gene names, but with bizzare ensembl ids, and shifted slightly in position. assuming these are bad.
    fullInfo = fullInfo[!(fullInfo$ensembl_gene_id %in% c("64245", "64244"))]
    return(fullInfo)
}


