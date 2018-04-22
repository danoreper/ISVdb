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


    sns = list()
    i = 1
    for(sn in c(115, 116))
    {
        
        sns[[i]] = fread(dat(paste0("b6_reference/mus_musculus.GRCm38.75/snord/",sn,".csv")))
        sns[[i]] = data.frame(ensembl_gene_id = as.character(sns[[i]]$id_with_url),
                              brainImprinted = F,
                              strainEffect = F,
                              Expressed.allele = NA)
        i = i + 1
    }
    
#collate the supplementary table imprinted genes, and the literature known imprinted genes into a single data frame
#we collect the bare minimum of what is necessary from these tables, and look up the rest from ensembl mart, as that is an online resource that is more current.
    impGenes =data.frame(
        ensembl_gene_id = c(as.character(litOnly$Ensembl.Gene.ID), as.character(imprinted$Ensembl.name)), 
        mgi_symbol      = c(as.character(litOnly$Gene), as.character(imprinted$Gene)),
        brainImprinted  = c(rep(F,nLit), rep(T, nImprint)),
        strainEffect    = c(rep(F,nLit), as.character(imprinted$Strain.effect.))=="Yes",
        Expressed.allele = c(as.character(litOnly$Expressed.allele), as.character(imprinted$Expressed.allele)))


    impGenes = data.table(impGenes)
    ## remove the bogus 115/116, lacking actual gene name, rows
    impGenes = impGenes[!mgi_symbol %in% c("Snord115","Snord116")]
    impGenes$mgi_symbol = NULL

    ##put the real snord records in
    for(sn in sns) {impGenes = rbind(impGenes, sn)}
    impGenes = data.table(impGenes,key="ensembl_gene_id")

    ## browser()

    
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
    
    ##some number of ensembl ids wont exist in ensembl mart, lets examine them
    allEnsembl   = unique(c(as.character(imprinted$Ensembl.name), as.character(literature$Ensembl.Gene.ID)))
    missingEnsembl = setdiff(allEnsembl, fullInfo$ensembl_gene_id)
    print("following ensembl ids were not found in biomart... they are deprecated as of my last manual check against ensembl page")
    print(missingEnsembl)

    ## fullInfo = rbind(fullInfo, snords, fill = T)
    #redundant gene names, but with bizzare ensembl ids, and shifted slightly in position. assuming these are bad.
    fullInfo = fullInfo[!(fullInfo$ensembl_gene_id %in% c("64245", "64244"))]
    return(fullInfo)
}


