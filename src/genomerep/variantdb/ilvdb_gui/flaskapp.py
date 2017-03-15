from flask import Flask, render_template, request, Response, stream_with_context
from model import InputForm, InputForm_gui
from query import generatequery, generatequery_gui
from urllib import quote,unquote

import MySQLdb
import MySQLdb.cursors # important

app = Flask(__name__)

def connection(): # Establish connection with ilvdb database/
    conn = MySQLdb.connect(host="valdardb.its.unc.edu",user = "valdar_user",passwd="u6GENwMGMPkDE6",db="exon_1410b",cursorclass = MySQLdb.cursors.SSCursor)
    c = conn.cursor()

    return c, conn

def generate_d(c): # The function that take a connection variable and streamming MyQTL to download file
    # input c. connection that already executed QTL query
    result_header= [x[0] for x in c.description]
    yield ','.join(str(e) for e in result_header) + '\n' # yield header
    
    row = c.fetchone()
    while row is not None:
        yield ','.join(str(e) for e in row) + '\n'  # yield data in a memoryless way
        row = c.fetchone()

@app.route('/large.csv')   # Route for testing download speed.
def generate_large_csv():  # Download may be interrupted by UNC firewall.
    c, conn = connection()
    c.execute("SELECT variant_id,chrom,pos,strain_name,allele_1,allele_2 from genotype_view limit 1000000")
    return Response(stream_with_context(generate_d(c)), mimetype='text/csv')

########################### The Main GUI System #########################
# Testing module 2
@app.route('/download/<name>') # Download route with <name> as the quoted version of query name
def download(name):
    c, conn = connection()
    try:
        c.execute(unquote(name))
        return Response(stream_with_context(generate_d(c)), mimetype='text/csv')
    except:
        return "Unknown query"

# Stable module
@app.route('/',methods=['GET','POST'])
def gui_page_test():
    c, conn = connection()
    form = InputForm_gui(request.form)
    if request.method == 'POST':
        result=generatequery_gui(variantinputtype = form.variantinputtype.data, view = form.view.data,
                                 variant_id =form.variant_id_text.data, variantfile = form.variantfile.data,
                                 chr = form.chr.data, strnum = form.strnum.data, endnum = form.endnum.data, bedfile = form.bedfile.data,
                                 genesearch=form.genesearch.data, genefile= form.genefile.data,
                                 selectCollum1=form.selectCollum1.data, selectCollum2=form.selectCollum2.data, selectCollum3=form.selectCollum3.data, selectCollum4= form.selectCollum4.data,
                                 prob_cutoff= form.prob_cutoff.data, ismax= form.ismax.data, homo = form.homo.data, consequence_list=form.consequence.data,
                                 strainlistnew = form.strainlistnew.data, strainp1new=form.strainp1new.data, strainp2new=form.strainp2new.data,
                                 limnum=form.limnum.data)
        form.variantinputtype.default=form.variantinputtype.data
        form.view.default = form.view.data
        c.execute(str(result))
        result_table = c.fetchall()
        result_header= [x[0] for x in c.description]
    else:
        result = "SELECT variant_id,chrom,pos,strain_name,allele_1,allele_2 from genotype_view limit 10"
        result_table = [[1L, '1', 3206310L, 'CC040_TauUncb38V01', 'AT', 'AT'], [1L, '1', 3206310L, 'CC001_Uncb38V01', 'AT', 'AT'], [1L, '1', 3206310L, 'CC002_Uncb38V01', 'AT', 'AT'], [1L, '1', 3206310L, 'CC003_Uncb38V01', 'AT', 'AT'], [1L, '1', 3206310L, 'CC009_Uncb38V01', 'AT', 'AT'], [1L, '1', 3206310L, 'CC012_GeniUncb38V01', 'AT', 'AT'], [1L, '1', 3206310L, 'CC013_GeniUncb38V01', 'AT', 'AT'], [1L, '1', 3206310L, 'CC026_GeniUncb38V01', 'AT', 'AT'], [1L, '1', 3206310L, 'CC037_TauUncb38V01', 'AT', 'AT'], [1L, '1', 3206310L, 'CC059_TauUncb38V01', 'AT', 'AT']]
        result_header = ['variant_id','chrom','pos','strain_name','allele_1','allele_2']
    result0 = result.split(" limit ",1)[0]  # The result0 query doesnot set limit
    result0q= quote(result0) # URLable version of query

    return render_template("index.html", form=form, result=result, result0q=result0q,result_table=result_table, result_header=result_header)


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

