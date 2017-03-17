import os, time, glob
from data import geneidchoice   #, genenamechoice

gene_id_t = {}
for item in geneidchoice:
    gene_id_t[item[1]] = item[0]


def generatequery(queryinput="SELECT distinct variant_id,chrom,pos,strain_name,allele_1,allele_2,prob,is_max,consequence_1,consequence_2,gene_name,transcript_name FROM genotype_transcript_view WHERE (strain_id=1 OR strain_id=2 OR strain_id=3 OR strain_id=4 OR strain_id=5 OR strain_id=6 OR strain_id=7 OR strain_id=8)"):
    if "limit" in queryinput:
	    return str(queryinput)
    else:
        return str(queryinput) + " limit 100 "


def generatequery_gui(variantinputtype='2', view='1',
                      variant_id="1,2,3",variantfile= None,
                      chr='ALL',strnum=0,endnum=3000000, bedfile=None,
                      genesearch = "", genefile= None,
                      selectCollum1 = [],selectCollum2 =[],selectCollum3= [],selectCollum4= [],
                      prob_cutoff = 0, ismax = 0,homo = 'ALL',consequence_list = [],
                      strainlistnew = [], strainp1new = [], strainp2new = [],
                      limnum = 100):
    strnum = int(strnum)
    endnum = int(endnum)
    
    FROM = ""
    SELECT = []
    
    consequence_list_b = []
    for x in consequence_list:
        for y in x.split(','):
            consequence_list_b.append(y)
    consequence_list = list(set(consequence_list_b))
    
    if variantinputtype=='1':
        if variantfile != 1:
            variant_list = variant_id.replace(" ", "").split(',')
        else:
            pass
    if variantinputtype=='2':
        if bedfile == None:
            [chr,strnum,endnum]=[chr,strnum,endnum]
    if variantinputtype == '3':
        if genefile != 1:
            gene_list = genesearch.replace(" ", "").split(',')
            gene_list = [gene_id_t[s] for s in gene_list]

    if view == '1':
        WHERE_list=[]
        selectCollum = selectCollum1
        FROM = 'exon_1410b.genotype_transcript_view'

        if 'consequence' in selectCollum:
            index = selectCollum.index('consequence')
            selectCollum = selectCollum[0:index]+ ["consequence_1","consequence_2"]+ selectCollum[index+1:]
        if 'allele' in selectCollum:
            index = selectCollum.index('allele')
            selectCollum = selectCollum[0:index]+ ["allele_1","allele_2"]+ selectCollum[index+1:]

        SELECT = ",".join(selectCollum)

        if homo == 'Homo':
            WHERE_list += ["allele_id_1 <=> allele_id_2"]
        if homo == 'Hetero':
            WHERE_list += ["NOT allele_id_1 <=> allele_id_2"]
        if len(consequence_list) > 0 :
            WHERE_list += ["( " + " OR ".join(["consequence_id_1=" + str(s) + " or consequence_id_2=" + str(s) for s in consequence_list]) + " )"]

        if variantinputtype=='1':
            WHERE_list += ["( " + " OR ".join(["variant_id = " + str(s) for s in variant_list]) + " )"]
        if variantinputtype=='2':
            if chr != "ALL":
                WHERE_list += ["chrom='"+str(chr)+"'"]
            if strnum != None:
                WHERE_list += ["pos>="+ str(strnum)]
            if endnum != None:
                WHERE_list += ["pos<="+ str(endnum)]
        if variantinputtype == '3':
            WHERE_list += ["( " + " OR ".join(["gene_id = " + str(s) for s in gene_list]) + " )"]

        if prob_cutoff>0:
            WHERE_list += ["prob>" + str(prob_cutoff)]
        if ismax:
            WHERE_list += ["is_max=1"]
    
        if len(strainlistnew)>0:
            WHERE_list += ["( " + " OR ".join(["strain_id = " + str(s) for s in strainlistnew]) + " )"]

        if len(WHERE_list)>0:
            WHERE = " and ".join(WHERE_list)
            query = "SELECT distinct "+ SELECT + " FROM " + FROM + " WHERE " + WHERE
        else:
            query = "SELECT distinct "+ SELECT + " FROM " + FROM

    ###################################################
    if view == '2':
        WHERE_list=[]
        selectCollum = selectCollum2
        FROM = 'exon_1410b.diplotype_view'

        if 'founder_name' in selectCollum:
            index = selectCollum.index('founder_name')
            selectCollum = selectCollum[0:index]+ ["founder_name_1","founder_name_2"]+ selectCollum[index+1:]
        
        SELECT = ",".join(selectCollum)

        if homo == 'Homo':
            WHERE_list += ["founder_name_1 <=> founder_name_2"]
        if homo == 'Hetero':
            WHERE_list += ["NOT founder_name_1 <=> founder_name_2"]

        if variantinputtype=='1':
            WHERE_list += ["( " + " OR ".join(["variant_id = " + str(s) for s in variant_list]) + " )"]
        if variantinputtype=='2':
            if chr != "ALL":
                WHERE_list += ["chrom='"+str(chr)+"'"]
            if strnum != None:
                WHERE_list += ["pos>="+ str(strnum)]
            if endnum != None:
                WHERE_list += ["pos<="+ str(endnum)]
        if variantinputtype == '3':
            WHERE_list += ["( " + " OR ".join(["gene_id = " + str(s) for s in gene_list]) + " )"]

        if prob_cutoff>0:
            WHERE_list += ["prob>" + str(prob_cutoff)]

        if len(strainlistnew)>0:
            WHERE_list += ["( " + " OR ".join(["strain_id = " + str(s) for s in strainlistnew]) + " )"]

        if len(WHERE_list)>0:
            WHERE = " and ".join(WHERE_list)
            query = "SELECT distinct "+ SELECT + " FROM " + FROM + " WHERE " + WHERE
        else:
            query = "SELECT distinct "+ SELECT + " FROM " + FROM


    ###################################################
    if view == '3':
        WHERE_list=[]
        selectCollum = selectCollum3
        FROM = "exon_1410b.genotype_sampling_view as s1 inner join exon_1410b.genotype_sampling_view as s2 on s1.variant_id=s2.variant_id and s1.transcript_id=s2.transcript_id";
        
        strain_pair = [(strainp1new,strainp2new)]

        SELECT = ",".join(selectCollum)

        if homo == 'Homo':
            WHERE_list += ["s1.allele_index <=> s2.allele_index"]
        if homo == 'Hetero':
            WHERE_list += ["NOT s1.allele_index <=> s2.allele_index"]

        if variantinputtype=='1':
            WHERE_list += ["( " + " OR ".join(["s1.variant_id = " + str(s) for s in variant_list]) + " )"]
            WHERE_list += ["( " + " OR ".join(["s2.variant_id = " + str(s) for s in variant_list]) + " )"]
        if variantinputtype=='2':
            if chr != "ALL":
                WHERE_list += ["s1.chrom='"+str(chr)+"'"]
                WHERE_list += ["s2.chrom='"+str(chr)+"'"]
            if strnum != None:
                WHERE_list += ["s1.pos>="+ str(strnum)]
                WHERE_list += ["s2.pos>="+ str(strnum)]
            if endnum != None:
                WHERE_list += ["s1.pos<="+ str(endnum)]
                WHERE_list += ["s2.pos<="+ str(endnum)]
        if variantinputtype == '3':
            WHERE_list += ["( " + " OR ".join(["s1.gene_id = " + str(s) for s in gene_list]) + " )"]
            WHERE_list += ["( " + " OR ".join(["s2.gene_id = " + str(s) for s in gene_list]) + " )"]

        if prob_cutoff>0:
            WHERE_list += ["s1.prob*s2.prob >" + str(prob_cutoff)]
        
        WHERE_list_strain = []
        for s2 in strain_pair:
            WHERE_list_strain +=  ["(s1.strain_id=" + str(s2[0]) + " AND s2.strain_id=" + str(s2[1]) + " )"]
        if len(consequence_list) > 0 :
            WHERE_list += ["( " + " OR ".join(["s1.consequence_id=" + str(s) + " OR s2.consequence_id=" + str(s) for s in consequence_list]) + " )"]

        WHERE_list += ["( " + " OR ".join([ s for s in WHERE_list_strain]) + " )"]


        if len(WHERE_list)>0:
            WHERE = " and ".join(WHERE_list)
            query = "SELECT "+ SELECT + " FROM " + FROM + " WHERE " + WHERE
        else:
            query = "SELECT "+ SELECT + " FROM " + FROM

    ###################################################
    if view == '4':
        WHERE_list=[]
        selectCollum = selectCollum4
        FROM = "exon_1410b.diplotype_sampling_view as s1 inner join exon_1410b.diplotype_sampling_view as s2 on s1.variant_id=s2.variant_id";
        
        strain_pair = [(strainp1new,strainp2new)]

        SELECT = ",".join(selectCollum)

        if homo == 'Homo':
            WHERE_list += ["s1.founder_id <=> s2.founder_id"]
        if homo == 'Hetero':
            WHERE_list += ["NOT s1.founder_id <=> s2.founder_id"]

        if variantinputtype=='1':
            WHERE_list += ["( " + " OR ".join(["s1.variant_id = " + str(s) for s in variant_list]) + " )"]
            WHERE_list += ["( " + " OR ".join(["s2.variant_id = " + str(s) for s in variant_list]) + " )"]
        if variantinputtype=='2':
            if chr != "ALL":
                WHERE_list += ["s1.chrom='"+str(chr)+"'"]
                WHERE_list += ["s2.chrom='"+str(chr)+"'"]
            if strnum != None:
                WHERE_list += ["s1.pos>="+ str(strnum)]
                WHERE_list += ["s2.pos>="+ str(strnum)]
            if endnum != None:
                WHERE_list += ["s1.pos<="+ str(endnum)]
                WHERE_list += ["s2.pos<="+ str(endnum)]
        if variantinputtype == '3':
            WHERE_list += ["( " + " OR ".join(["s1.gene_id = " + str(s) for s in gene_list]) + " )"]
            WHERE_list += ["( " + " OR ".join(["s2.gene_id = " + str(s) for s in gene_list]) + " )"]

        if prob_cutoff>0:
            WHERE_list += ["s1.prob*s2.prob >" + str(prob_cutoff)]
        
        WHERE_list_strain = []
        for s2 in strain_pair:
            WHERE_list_strain +=  ["(s1.strain_id=" + str(s2[0]) + " AND s2.strain_id=" + str(s2[1]) + " )"]

        WHERE_list += ["( " + " OR ".join([ s for s in WHERE_list_strain]) + " )"]


        if len(WHERE_list)>0:
            WHERE = " and ".join(WHERE_list)
            query = "SELECT "+ SELECT + " FROM " + FROM + " WHERE " + WHERE
        else:
            query = "SELECT "+ SELECT + " FROM " + FROM
    if (int(limnum)==0):
        return str(query)
    else:
        return str(query) + " limit " + str(int(limnum))

if __name__ == '__main__':
    print generatequery_gui(variantinputtype='1', view='1', variant_id="1",variantfile= None, chr='ALL',strnum=30000000,endnum=31000000, bedfile=None, genesearch = "", genefile= None, selectCollum1 =["variant_id","chrom","pos","strain_name","allele","prob","is_max","consequence","gene_name","transcript_name"],selectCollum2 =[],selectCollum3= [],selectCollum4= [], prob_cutoff = 0.5, ismax = 1,homo = 'ALL',consequence_list = ['1','2'], strainlistnew = ['1','2','3','4','5','6','7','8'], strainp1new = [], strainp2new = [], limnum = 100)

