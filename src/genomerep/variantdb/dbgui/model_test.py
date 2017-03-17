from wtforms import widgets,Form,validators,fields
from data import strain_list,consequence_list
consequence_list_s = [s[1] for s in consequence_list]
consequence_list = zip(consequence_list_s,consequence_list_s)

chrom_list_choice=zip(range(1,20),range(1,20))  # zip(range(1,20),range(1,20))
chrom_list_choice=chrom_list_choice +  [('X','X'),('Y','Y'),('MT','MT'),('ALL','ALL')] #+[

genotype_view_choice_seed=["variant_id","chrom","pos","strain_name","allele","prob","is_max","gene_name","transcript_name","consequence"]
genotype_view_choice=zip(genotype_view_choice_seed,genotype_view_choice_seed)

diplotype_view_choice_seed=["variant_id","chrom","pos","strain_name","founder_name","prob","gene_name"]
diplotype_view_choice=zip(diplotype_view_choice_seed,diplotype_view_choice_seed)

genotype_cross_view_seed1=["variant_id","chrom","pos","strain_name","allele","prob","gene_name", "transcript_name","consequence"]
genotype_cross_view_seed2=["s1.variant_id as variant_id","s1.chrom as chrom","s1.pos as pos", "s1.strain_name as strain_name_1, s2.strain_name as strain_name_2", "s1.allele as allele_1, s2.allele allele_2","s1.prob*s2.prob as prob","s1.gene_name as gene_name","s1.transcript_name as transcript_name","s1.consequence_name as consequence_1, s2.consequence_name as consequence_2"]
genotype_cross_view_choice=zip(genotype_cross_view_seed2,genotype_cross_view_seed1)

diplotype_cross_view_seed1=["variant_id","chrom","pos","strain_name","founder_name","prob","gene_name"]
diplotype_cross_view_seed2=["s1.variant_id as variant_id","s1.chrom as chrom","s1.pos as pos","s1.strain_name as strain_1, s2.strain_name as strain_2","s1.founder_name as haplotype_1, s2.founder_name as haplotype_2","s1.prob*s2.prob as prob", "s1.gene_name as gene_name"]
diplotype_cross_view_choice=zip(diplotype_cross_view_seed2,diplotype_cross_view_seed1)

viewchoices=[('1', 'Genotype'), ('2', 'Diplotype'), ('3', 'Genotype cross'),('4','Diplotype cross')]


class InputForm_gui_test(Form):
    # No Interactive selection version. 2016.07.28	
    variantinputtype = fields.SelectField(label='In ', choices=[('1', 'Variant IDs'), ('2', 'Region(bp)'), ('3', 'Genes')],default='2')
    
	########## Specify method of selection ########
    variant_id_text= fields.TextAreaField(label='variant_id_text',default='Comma separated ISVDB variant IDs')
    #variantfile = fields.FileField('Upload variant_id file(Unavailable)')

    chr = fields.SelectField(label='Chrom', choices=chrom_list_choice,default='11')
    strnum = fields.FloatField(label='Start', default=30000000)
    endnum = fields.FloatField(label='End', default=35000000)
    #bedfile = fields.FileField('Upload bed file (Unavailable)')
   
    genesearch = fields.TextAreaField(label='Comma separated Ensembl IDs',default='ENSMUSG00000051951, ENSMUSG00000025931, ENSMUSG00000026596')
    #gene_list = fields.SelectField(label='Genes search result', choices=genenamechoice)
    #genefile = fields.FileField('Upload gene list file(Unavailable)')
	######### Specify view table to use #########
    view = fields.SelectField(label='View ', choices=viewchoices)
    # currently just use all selection?
    selectCollum1 = fields.SelectMultipleField(label='Select included fields',choices=genotype_view_choice,default=genotype_view_choice_seed)
    selectCollum2 = fields.SelectMultipleField(label='Select included fields',choices=diplotype_view_choice,default=diplotype_view_choice_seed)
    selectCollum3 = fields.SelectMultipleField(label='Select included fields',choices=genotype_cross_view_choice,default=genotype_cross_view_seed2)
    selectCollum4 = fields.SelectMultipleField(label='Select included fields',choices=diplotype_cross_view_choice,default=diplotype_cross_view_seed2)

	######## Specify target strain/founder ########
    strainlistnew = fields.SelectMultipleField(label='Of strains',choices= strain_list, default=['3','8','9','76'])
    strainp1new   = fields.SelectField(label=' of S1',choices= strain_list, default='8')
    strainp2new   = fields.SelectField(label=' with S2',choices= strain_list, default='46')
    
	######## add filter to final data #########
    prob_cutoff = fields.FloatField(label='Prob >=', default=0.00)
    homo = fields.SelectField(label='Zygosity:', choices=[('Homo', 'Homozygous'), ('Hetero', 'Heterozygous'), ('All', 'All')],default='All')
    ismax = fields.BooleanField(label='Only show max probability?',default=0)
    consequence = fields.SelectMultipleField(label='With consequences',choices= consequence_list)

    ############### limit
    limnum = fields.FloatField(label='Output limit', default=1000)





