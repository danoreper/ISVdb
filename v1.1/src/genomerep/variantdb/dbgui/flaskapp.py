import copy

from flask import Flask, render_template, request, Response, stream_with_context
from model import InputForm, InputForm_gui
from query import generatequery, generatequery_gui
from urllib import quote,unquote

from model_test import InputForm_gui_test
from query_test import generatequery_gui_test

# These two modules exist in the server but not local
import MySQLdb
import MySQLdb.cursors # important

app = Flask(__name__)

def connection(): # Establish connection with ilvdb database/
    conn = MySQLdb.connect(host="valdardb.its.unc.edu",user = "valdar_user",db="exon_1410E_11",passwd="u6GENwMGMPkDE6",cursorclass = MySQLdb.cursors.SSCursor)
    c = conn.cursor()

    return c, conn

def tablename(s): # change names in the result_header
    if 'gene_name' in s:
        s[s.index('gene_name')]='gene_id'
    if 'transcript_name' in s:
        s[s.index('transcript_name')]='transcript_id'
    if 'founder_name_1' in s:
        s[s.index('founder_name_1')]='haplotype_1'
    if 'founder_name_2' in s:
        s[s.index('founder_name_2')]='haplotype_2'
    if 'strain_name' in s:
        s[s.index('strain_name')]='strain'
    if 'strain_name_1' in s:
        s[s.index('strain_name_1')]='strain_1'
    if 'strain_name_2' in s:
        s[s.index('strain_name_2')]='strain_2'
    return(s)

def generate_d(c): # The function that take a connection variable and streamming MyQTL to download file
    # input c. connection that already executed QTL query
    result_header= [x[0] for x in c.description]
    result_header=tablename(result_header)
    yield ','.join(str(e) for e in result_header) + '\n' # yield header
    
    row = c.fetchone()
    while row is not None:
        yield ','.join(str(e) for e in row) + '\n'  # yield data in a memoryless way
        row = c.fetchone()

def strain_sep(result_table, strain_id):
    for i in range(0,len(result_table)):
        a = list(result_table[i])
        strnld=str(a[strain_id])
        if strnld=='A_J' or strnld =='CAST_EiJ' or strnld =='129S1_SvlmJ' or strnld =='NOD_ShiLtJ' or strnld =='NZO_HlLtJ' or strnld =='PWK_PhJ' or strnld =='WSB_EiJ' or strnld =='C57BL6J':
            continue
        else:
            a[strain_id] = strnld.split('_')[0]
            result_table[i] = a
    return result_table

########################### The Main GUI System #########################
@app.route('/download/<name>') # Download route with <name> as the quoted version of query name
def download(name):
    c, conn = connection()
    try:
        c.execute(unquote(name))
        return Response(stream_with_context(generate_d(c)), mimetype='text/csv')
    except:
        return "Unknown query"

# Stable module
@app.route('/v1.0',methods=['GET','POST'])
def gui_page():
    c, conn = connection()
    form = InputForm_gui(request.form)
    if request.method == 'POST':
        result=generatequery_gui(variantinputtype = form.variantinputtype.data, view = form.view.data,
                                 variant_id =form.variant_id_text.data,
                                 chr = form.chr.data, strnum = form.strnum.data, endnum = form.endnum.data,
                                 genesearch=form.genesearch.data,
                                 selectCollum1=form.selectCollum1.data, selectCollum2=form.selectCollum2.data, selectCollum3=form.selectCollum3.data, selectCollum4= form.selectCollum4.data,
                                 prob_cutoff= form.prob_cutoff.data, ismax= form.ismax.data, homo = form.homo.data, consequence_list=form.consequence.data,
                                 strainlistnew = form.strainlistnew.data, strainp1new=form.strainp1new.data, strainp2new=form.strainp2new.data,
                                 limnum=form.limnum.data)
        form.variantinputtype.default=form.variantinputtype.data
        form.view.default = form.view.data
    else:
        result = "SELECT distinct variant_id,chrom,pos,strain_name,allele_1,allele_2,prob,is_max,consequence_1,consequence_2,gene_name,transcript_name FROM genotype_transcript_view WHERE chrom='19' and pos>=30000000 and pos<=35000000 and ( strain_id = 1 OR strain_id = 2 OR strain_id = 3 OR strain_id = 4 OR strain_id = 5 OR strain_id = 6 OR strain_id = 7 OR strain_id = 8 OR strain_id = 9 ) limit 100";
    c.execute(str(result))
    result_table = c.fetchall()
    result_header= [x[0] for x in c.description]
    result0 = result.split(" limit ",1)[0]  # The result0 query doesnot set limit
    result0q= quote(result0) # URLable version of query
    result_header=tablename(result_header)

    return render_template("index.html", form=form, result=result, result0q=result0q,result_table=result_table, result_header=result_header)

@app.route('/help/',methods=['GET','POST'])
def help_my():
    return render_template("help_ilvdb.html")


@app.route('/',methods=['GET','POST'])
def gui_page_test():
    
    form = InputForm_gui_test(request.form)
    c, conn = connection()
    if request.method == 'POST':
        result=generatequery_gui_test(variantinputtype = form.variantinputtype.data, view = form.view.data,
                                 variant_id =form.variant_id_text.data,
                                 chr = form.chr.data, strnum = form.strnum.data, endnum = form.endnum.data,
                                 genesearch=form.genesearch.data,
                                 selectCollum1=form.selectCollum1.data, selectCollum2=form.selectCollum2.data, selectCollum3=form.selectCollum3.data, selectCollum4= form.selectCollum4.data,
                                 prob_cutoff= form.prob_cutoff.data, ismax= form.ismax.data, homo = form.homo.data, consequence_list=form.consequence.data,
                                 strainlistnew = form.strainlistnew.data, strainp1new=form.strainp1new.data, strainp2new=form.strainp2new.data,
                                 limnum=form.limnum.data)
    else:
        result = "SELECT distinct variant_id,chrom,pos,strain_name,allele_1,allele_2,prob,is_max,consequence_1,consequence_2,gene_name,transcript_name FROM exon_1410D_11.genotype_transcript_view WHERE chrom='11' and pos>=30000000 and pos<=35000000 and ( strain_id = 1 OR strain_id = 2 OR strain_id = 3 OR strain_id = 4 OR strain_id = 5 OR strain_id = 6 OR strain_id = 7 OR strain_id = 8 OR strain_id = 9 ) limit 100";

    c.execute(str(result))
    result_table = list(c)
    result_header= [x[0] for x in c.description]

    result0 = result.split(" limit ",1)[0]  # The result0 query doesnot set limit
    result0q= quote(result0) # URLable version of query
    result_header=tablename(result_header)

    result_table_b = copy.deepcopy(result_table) # Important! Must copy to list tuple to change it!
    if 'prob' in result_header: # chop prob to 0.000 total 5 charactors
        prob_id = result_header.index('prob')
        for i in range(0,len(result_table)):
            a = list(result_table_b[i])
            a[prob_id] = str(a[prob_id])[:5]
            result_table_b[i]=a
    if 'strain' in result_header:
        result_table = strain_sep(result_table_b, result_header.index('strain'))
    if 'strain_1' in result_header:
        result_table = strain_sep(result_table_b, result_header.index('strain_1'))
    if 'strain_2' in result_header:
        result_table = strain_sep(result_table_b, result_header.index('strain_2'))


    return render_template("index_test.html", form=form, result=result, result0q=result0q, result_table=result_table_b, result_header=result_header)


# fossial
################################ The copy and paste GUI system for debug #################
@app.route('/query/',methods=['GET','POST'])
def index():
    c, conn = connection()
    form = InputForm(request.form)
    if request.method == 'POST':
        if form.validate():
            result=generatequery(form.query_string.data)
            c.execute(str(result))
            result_table = c.fetchall()
            result_header= [x[0] for x in c.description]
        else:
            result='Invalidate Form'
            result_table = None
            result_header= None
    else:
        result = "SELECT variant_id,chrom,pos,strain_name,allele_1,allele_2 from genotype_view limit 10"
        result_table = [[1L, '1', 3206310L, 'CC040_TauUncb38V01', 'AT', 'AT'], [1L, '1', 3206310L, 'CC001_Uncb38V01', 'AT', 'AT'], [1L, '1', 3206310L, 'CC002_Uncb38V01', 'AT', 'AT'], [1L, '1', 3206310L, 'CC003_Uncb38V01', 'AT', 'AT'], [1L, '1', 3206310L, 'CC009_Uncb38V01', 'AT', 'AT'], [1L, '1', 3206310L, 'CC012_GeniUncb38V01', 'AT', 'AT'], [1L, '1', 3206310L, 'CC013_GeniUncb38V01', 'AT', 'AT'], [1L, '1', 3206310L, 'CC026_GeniUncb38V01', 'AT', 'AT'], [1L, '1', 3206310L, 'CC037_TauUncb38V01', 'AT', 'AT'], [1L, '1', 3206310L, 'CC059_TauUncb38V01', 'AT', 'AT']]
        result_header = ['variant_id','chrom','pos','strain_name','allele_1','allele_2']
    return render_template("index2.html", form=form, result=result, result_table=result_table, result_header=result_header)

if __name__ == "__main__":
	app.run(debug=True)

