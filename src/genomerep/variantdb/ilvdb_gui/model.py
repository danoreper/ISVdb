from wtforms import widgets,Form,validators,fields
from data import strain_list,consequence_list

chrom_list_choice=zip(range(1,20),range(1,20))
chrom_list_choice=chrom_list_choice+[('X','X'),('ALL','ALL')]

genotype_view_choice_seed=["variant_id","chrom","pos","strain_name","allele","prob","is_max","consequence","gene_name","transcript_name"]
genotype_view_choice=zip(genotype_view_choice_seed,genotype_view_choice_seed)

diplotype_view_choice_seed=["variant_id","chrom","pos","founder_name","prob","gene_name"]
diplotype_view_choice=zip(diplotype_view_choice_seed,diplotype_view_choice_seed)

genotype_cross_view_seed1=["variant_id","chrom","pos","gene_name","strain_name","allele","prob","transcript_name","consequence"]
genotype_cross_view_seed2=["s1.variant_id as variant_id","s1.chrom as chrom","s1.pos as pos","s1.gene_name as gene_name", "s1.strain_name as strain_name_s1, s2.strain_name as strain_name_s2", "s1.allele as allele_s1, s2.allele allele_s2","s1.prob*s2.prob as prob","s1.transcript_name as transcript_name","s1.consequence_name as consequence_s1, s2.consequence_name as consequence_s2"]
genotype_cross_view_choice=zip(genotype_cross_view_seed2,genotype_cross_view_seed1)

diplotype_cross_view_seed1=["variant_id","chrom","pos","gene_name","strain_name","founder_name","prob"]
diplotype_cross_view_seed2=["s1.variant_id as variant_id","s1.chrom as chrom","s1.pos as pos","s1.gene_name as gene_name","s1.strain_name as strain_name_s1, s2.strain_name as strain_name_s2","s1.founder_name as founder_s1, s2.founder_name as founder_s2","s1.prob*s2.prob as prob"]
diplotype_cross_view_choice=zip(diplotype_cross_view_seed2,diplotype_cross_view_seed1)

viewchoices=[('1', 'Genotype view'), ('2', 'Diplotype view'), ('3', 'Genotype cross view'),('4','Diplotype cross view')]

class InputForm(Form):
    query_string= fields.TextAreaField(label='Input Query',validators=[validators.InputRequired()])

class InputForm_gui(Form):
    # No Interactive selection version. 2016.07.28	
    variantinputtype = fields.SelectField(label='Select input type', choices=[('1', 'Variant_id'), ('2', 'Region'), ('3', 'Gene')],default='2')
    
	########## Specify method of selection ########
    variant_id_text= fields.TextAreaField(label='Input specific variant_id (Comma seperated)')
    variantfile = fields.FileField('Upload variant_id file(Unavailable)')

    chr = fields.SelectField(label='Chromosome', choices=chrom_list_choice,default='19')
    strnum = fields.FloatField(label='Start site (bp)', default=30000000)
    endnum = fields.FloatField(label='End site (bp)', default=35000000)
    bedfile = fields.FileField('Upload bed file (Unavailable)')
   
    genesearch = fields.TextAreaField(label='Gene Ensembl ID (Comma seperated)')
    #gene_list = fields.SelectField(label='Genes search result', choices=genenamechoice)
    genefile = fields.FileField('Upload gene list file(Unavailable)')
	######### Specify view table to use #########
    view = fields.SelectField(label='Select view type', choices=viewchoices)
    # currently just use all selection?
    selectCollum1 = fields.SelectMultipleField(label='Select included fields',choices=genotype_view_choice,default=genotype_view_choice_seed)
    selectCollum2 = fields.SelectMultipleField(label='Select included fields',choices=diplotype_view_choice,default=diplotype_view_choice_seed)
    selectCollum3 = fields.SelectMultipleField(label='Select included fields',choices=genotype_cross_view_choice,default=genotype_cross_view_seed2)
    selectCollum4 = fields.SelectMultipleField(label='Select included fields',choices=diplotype_cross_view_choice,default=diplotype_cross_view_seed2)

	######## Specify target strain/founder ########
    strainlistnew = fields.SelectMultipleField(label='Strains/Founders',choices= strain_list, default=['1','2','3','4','5','6','7','8'])
    strainp1new   = fields.SelectField(label='Strain1 in cross',choices= strain_list)
    strainp2new   = fields.SelectField(label='Strain2 in cross',choices= strain_list)
    
	######## add filter to final data #########
    prob_cutoff = fields.FloatField(label='Prob >=', default=0.00)
    homo = fields.SelectField(label='Heterozygous?', choices=[('Homo', 'Homozygous'), ('Hetero', 'Heterozygous'), ('All', 'All')],default='All')
    ismax = fields.BooleanField(label='Probablity is max?',default=1)
    consequence = fields.SelectMultipleField(label='Consequence? (default is All))',choices= consequence_list)

    ############### limit
    limnum = fields.FloatField(label='Output limit', default=100)





