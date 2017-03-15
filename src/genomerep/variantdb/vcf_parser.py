DEBUG = False

import os
# if(DEBUG):
#     os.chdir("../")
import gzip; 
import time;

import argparse;
# import sqlite3;
import MySQLdb as mdb;
from collections import OrderedDict
import genomerep.variantdb.postProcessParsedDB
import logging.config
import yaml
import ipdb


logging.config.dictConfig(yaml.load(open("../config/logging.config",'r')))
logger=logging.getLogger("mylog")
logger.debug("logging!")

    
INSERT_BATCH_SIZE = 5000
##engine   = "MYISAM"  

def getInputColIndsForOutput(headerTokens, outputFieldToInputField):
    colind = OrderedDict()
    for tableField in outputFieldToInputField.keys():
        headerField = outputFieldToInputField[tableField]
        ind = headerTokens.index(headerField)
        colind[tableField] = ind
    return colind

def navigateToheaderTokens(fle):
    logger.debug("navigating to header tokens")
    for line in fle:
        if (not line.startswith("##")):
            headerTokens = line.strip().split("\t")
            break

    logger.debug("header tokens are: " + str(headerTokens))
    return headerTokens


def getFirstHeader(vcfFiles):
    logger.debug("trying to pull out header tokens")
    logger.debug("current path is: " + os.getcwd())
    logger.debug("opening:" + vcfFiles[0])
    ##grab the first vcf, we assume all vcfs will have the same header tokens. Will check this later on
    fle = gzip.open(vcfFiles[0], 'r')
    logger.debug("opened file")
    headerTokensGlobal = navigateToheaderTokens(fle)
    fle.close()
    return headerTokensGlobal

def dropAllTables(con):
    con.query("Drop table if exists inbred_line_founder;")
    con.query("Drop table if exists consequence_info;")
    con.query("Drop table if exists founder_genotype;")
    con.query("Drop table if exists genotype;")
    con.query("Drop table if exists genotype_sampling;")
    con.query("Drop table if exists inbred_line_fixed;")
    con.query("Drop table if exists strain;")
    con.query("Drop table if exists consequence;")
    con.query("Drop table if exists allele;")
    con.query("Drop table if exists gene;")
    con.query("Drop table if exists transcript;")
    con.query("Drop table if exists info;")
    con.query("Drop table if exists variant;")

def parsevcf(vcfFiles, foundersFile, dbname, host, user, passwd, firstVariantId, allowedTranscriptsFile, limit):
    logger.debug("!!!!!!!!!!!!!!!!!!!!!")
    logger.debug("about to begin parsing for db " + dbname)
    tmStart = time.time()
    con = mdb.connect(host=host, user=user, passwd=passwd)
    cursor = con.cursor()
    createcommand ="CREATE DATABASE IF NOT EXISTS "+dbname +";"
    logger.debug(createcommand) 
    cursor.execute(createcommand)
    cursor.execute("USE "+dbname+";")
    con.query("set global max_allowed_packet=8048576000;")

    allowedTranscripts = []
    if not allowedTranscriptsFile is None:
        fh = open(allowedTranscriptsFile)
        for transcript in fh:
            allowedTranscripts.append(transcript.strip())

   #ipdb.set_trace()           
    if limit is None:
        limit = float("inf")

    logger.debug("limit is" + str(limit))
    flag = True    
    if flag:
        dropAllTables(con)
        logger.debug("dropped all tables")
        import glob;
        vcfFiles = os.path.join(vcfFiles,"*vcf.gz")
        logger.debug("reading from " + vcfFiles)
        vcfFiles = glob.glob(vcfFiles)
        logger.debug(vcfFiles)
        headerTokensGlobal = getFirstHeader(vcfFiles)
        logger.debug("got header tokens: ") 
        
        diploParser = DiploParsing(headerTokensGlobal, foundersFile)
        logger.debug("got diplo parser")
        parsers = [ 
                    VariantParsing(headerTokensGlobal) ,
                    AlleleParsing(headerTokensGlobal),
                    diploParser,
                    DiploMetadataParsing(headerTokensGlobal, foundersFile),
                    InfoParsing(headerTokensGlobal)
                  ] 
        logger.debug("got all other parsers")
    
        writers = [parser.getDbWriter(con) for parser in parsers]    
        logger.debug("got all other writers")

        variant_id = 1
        if(not firstVariantId is None):
            variant_id= int(firstVariantId)
            
        logger.debug("about to iterate over vcffiles")
        for vcfFile in vcfFiles:
            logger.debug("parsing " + vcfFile) 
            if variant_id>limit:
                break
            logger.debug(vcfFile)
            vcfFile = gzip.open(vcfFile, 'r')
            headerTokens = navigateToheaderTokens(vcfFile)
            if ",".join(headerTokens)!=",".join(headerTokensGlobal):
                raise("inconsistent headers between vcf files")
            
            for line in vcfFile:

                if variant_id>limit:
                    break
                ##logger.debug(line)
                vcfTokens = line.split("\t")
                if diploParser.isValidVariant(vcfTokens):
                    for i in range(len(parsers)):
                        parser = parsers[i]
                        writer = writers[i]
                        values = parser.getValues(vcfTokens, variant_id)
                        if not values==None:
                            if DEBUG:
                                try:
                                    writer.write(values)
                                except:
                                    ipdb.set_trace()
                                    print(5)
                            else:
                                writer.write(values)
                variant_id = variant_id + 1
                if (variant_id % 10000)==0:
                    print(variant_id)
            vcfFile.close()

        #close writers only after we are done with all vcfs
        [writerr.close() for writerr in writers]
    
    logger.debug("attempting to postprocess")
    genomerep.variantdb.postProcessParsedDB.postprocess(con, allowedTranscripts, engine)
    tmEnd = time.time();
    logger.debug(tmEnd-tmStart)
    con.close()
class VariantParsing:
    
    _outputFieldToInputField = OrderedDict([("chrom" , "#CHROM"),
                                            ("pos"   , "POS")])  
    _outputfields = ["variant_id"]
    _outputfields.extend(_outputFieldToInputField.keys())
     

    def __init__(self, headerTokens):      
        self.colind = getInputColIndsForOutput(headerTokens, VariantParsing._outputFieldToInputField).values()
        
    def getValues(self, vcftokens, variant_id):
        value = [variant_id];
        value.extend([vcftokens[ind] for ind in self.colind])
        values = [value] #2d list of lists, even though theres only a single set of values.
        return(values)
        
    @staticmethod
    def getDbWriter(con):
        tablename = "variant"
        VariantParsing.buildTable(con, tablename)

        writer = TableWriter(con, VariantParsing._outputfields, tablename, INSERT_BATCH_SIZE)
        return(writer)

    @staticmethod
    def buildTable(con, tablename):
#         con.query("DROP TABLE IF EXISTS "+ tablename + ";")
        con.query('''CREATE TABLE IF NOT EXISTS ''' + tablename + '''
                   (variant_id integer unsigned AUTO_INCREMENT PRIMARY KEY NOT NULL, chrom varchar(6) NOT NULL, pos integer NOT NULL) ENGINE = ''' + engine + ''';''')


class DiploParsing:
    
    @staticmethod
    def _parseFounders(foundersFile):
        founders = []
        founder_file = open(foundersFile, 'rb')
        for row in founder_file:
            fndr = row.strip()
            if(fndr!="C57BL6J"):
                founders.append(row.strip())
        founder_file.close()
        return founders    
        
    _outputfields = ["variant_id","strain_name","allele_index_1","allele_index_2","prob", "is_max" ]   
#     _outputfields = ["variant_id","strain_name","allele_index", "prob"]
    def __init__(self, headerTokens, foundersFile):

        #see vcf spec for gl: we are assuming biallelic genotypes of form j/k, and no more than 200 alleles...
        maxAlleleInd       = 200
        numConditions      = (maxAlleleInd*(maxAlleleInd+1)/2)+maxAlleleInd
        self.jalleles_by_F = numConditions*[-1]
        self.kalleles_by_F = numConditions*[-1]
        secondaryTokens    = "GT:GQ:DP:MQ0F:GP:PL:AN:MQ:DV:DP4:SP:SGB:PV4:FI"
        secondaryTokens    = secondaryTokens.split(":")
        self.gp_ind        = secondaryTokens.index("GP")

        for j in range(0,maxAlleleInd-1):
            for k in range(j,maxAlleleInd):
                F = (k*(k+1)/2)+j 
                self.jalleles_by_F[F] = j
                self.kalleles_by_F[F] = k
            
        founders  = DiploParsing._parseFounders(foundersFile)
        
        self.founderToInt = {}
        counter = 1
        for founder in founders:
            if not founder in self.founderToInt.keys():
                self.founderToInt[founder]=counter
                counter = counter + 1 
            
        #building this map just to use the convenient method for finding header column indexes
        outputFieldToInputField = OrderedDict()
        for founder in founders:
            outputFieldToInputField[founder] = founder
        self.founderinds = getInputColIndsForOutput(headerTokens, outputFieldToInputField)
    
    def isValidVariant(self, vcftokens):
        for colind in self.founderinds.itervalues():
            diplotypeVal = vcftokens[colind].split(":")
            diplotypeVal = diplotypeVal[0]
            if diplotypeVal!="0/0" and diplotypeVal!=".":
                return True
        else: 
            return False;    
            

    def getProbs(self, gps):
        probs = len(gps) * [0]
        for i in range(0, len(gps)):
            gp = gps[i]
            if(gp=="."):
                probs[i] = 0
            else:
                probs[i] = pow(10, int(gp) / -10)
        return probs

    def getValues(self, vcftokens, variant_id):
        values = []
        for founder, colind in self.founderinds.iteritems():
            if(len(vcftokens)<=colind):
                logger.debug("warn!")
                logger.debug(vcftokens)
                return(None)
            diplotypeVal = vcftokens[colind].split(":")
            
            gps = diplotypeVal[self.gp_ind].split(",")
            mleDiplotype = diplotypeVal[0] 
#             kk = diplotypeVal[0]
#             kk = kk.split("/")
#             if (kk[0]=="." and kk[1]!=".") or (kk[1]=="." and kk[0]!="."):
#                 logger.debug("err")
            
            if mleDiplotype=="./.":
                values.append([variant_id, founder, None, None, 1, True])
            elif all([x == "." for x in gps]):
                values.append([variant_id, founder, 0, 0, 1, True])
            else:
                probs = self.getProbs(gps)
                probSum = float(sum(probs))
                probs = [prob/probSum for prob in probs]
                maxProbIndex = probs.index(max(probs))
                
                for i in range(0,len(gps)):
                    if probs[i]>=.01:
                        allele_index_1 = self.jalleles_by_F[i]
                        allele_index_2 = self.kalleles_by_F[i]
                        values.append([variant_id, 
                                       founder, 
                                       allele_index_1, 
                                       allele_index_2, 
                                       probs[i], 
                                       i==maxProbIndex])
                    
        return(values)
            
    @classmethod
    def getDbWriter(cls,con):
        #building a diplotype writer        
        tablename = "founder_genotype"
        cls.buildTable(con, tablename)
        writer = TableWriter(con, cls._outputfields, tablename, INSERT_BATCH_SIZE)
        return(writer)

    @staticmethod
    def buildTable(con, tablename):
        #TODO maxchar
        maxchar = 30
        con.query("DROP TABLE IF EXISTS " +tablename +";")
        con.query('''CREATE TABLE IF NOT EXISTS ''' + tablename + '''
                    (variant_id integer unsigned NOT NULL, 
                    strain_name varchar('''+str(maxchar)+''') NOT NULL, 
                    allele_index_1 tinyint unsigned, 
                    allele_index_2 tinyint unsigned, 
                    prob float unsigned NOT NULL,
                    is_max tinyint unsigned NOT NULL,
                    
                    FOREIGN KEY(variant_id) REFERENCES variant(variant_id),
                    FOREIGN KEY(allele_index_1)   REFERENCES allele(allele_index),
                    FOREIGN KEY(allele_index_2)   REFERENCES allele(allele_index)
                    ) ENGINE = ''' + engine + ''';''')



class DiploMetadataParsing:
        
    _outputInfoFieldToInputField = OrderedDict([("gq" ,  "GQ"),
                                            ("dp" ,  "DP"),
                                            ("mq0f", "MQ0F"),
                                            ("an",   "AN"),
                                            ("mq",   "MQ"),
                                            ("dv",   "DV"),
                                            ("dp4",  "DP4"),
                                            ("sp",   "SP"),
                                            ("sgb",  "SGB"),
#                                             ("pv4",  "PV4"),
                                            ("fi",    "FI")])  
    
    #gt and gp not included here because they are already taken into account in the DiploParsing
    _outputfields = ["variant_id","strain_name"]
    _outputfields.extend(_outputInfoFieldToInputField.keys())
    
    insertInd     = _outputfields.index("dp4")
    _outputfields.remove("dp4")
    _outputfields.insert(insertInd,  "dp_ref_fwd")
    _outputfields.insert(insertInd+1,"dp_ref_rev")
    _outputfields.insert(insertInd+2,"dp_alt_fwd")
    _outputfields.insert(insertInd+3,"dp_alt_rev")

    secondaryTokens    = "GT:GQ:DP:MQ0F:GP:PL:AN:MQ:DV:DP4:SP:SGB:PV4:FI"
#         "GQ: DP: MQ0F: GP:       AN: MQ: DV: DP4:     SP: SGB:       PV4: FI"
#         15:  4:  0.25: 80,15,0:  2:  25: 4:  0,0,4,0: 0:  -0.556411: .:   0

    secondaryTokens       = secondaryTokens.split(":")
    _infoinds         = getInputColIndsForOutput(secondaryTokens, _outputInfoFieldToInputField)
    
    def __init__(self, headerTokens, foundersFile):
        #see vcf spec for gl: we are assuming biallelic genotypes of form j/k, and no more than 200 alleles...
        #maxAlleleInd       = 200
        founders  = DiploParsing._parseFounders(foundersFile)
        #building this map just to use the convenient method for finding header column indexes
        outputFieldToInputField = OrderedDict()
        for founder in founders:
            outputFieldToInputField[founder] = founder
        self.founderinds = getInputColIndsForOutput(headerTokens, outputFieldToInputField)
        
    def getValues(self, vcftokens, variant_id):
        values = []
        for founder, colind in self.founderinds.iteritems():
            if(len(vcftokens)<=colind):
                logger.info("warn!")
                logger.info(vcftokens)
                return(None)
            
            infotokens = vcftokens[colind].split(":")
            infotokens = [tok.strip() for tok in infotokens]
            _infoinds = DiploMetadataParsing._infoinds
            dp4 = infotokens[_infoinds["dp4"]]
            dp4 = dp4.split(",")
            if(len(dp4)==1):
                dp4 = [".", ".", ".", "."] 
            
            value = [variant_id, founder] 
            valueExtra = (infotokens[_infoinds["gq"]],
                    infotokens[_infoinds["dp"]],
                    infotokens[_infoinds["mq0f"]],
                    infotokens[_infoinds["an"]],
                    infotokens[_infoinds["mq"]],
                    infotokens[_infoinds["dv"]],
                    dp4[0],
                    dp4[1],
                    dp4[2],
                    dp4[3],
                    infotokens[_infoinds["sp"]],
                    infotokens[_infoinds["sgb"]],
                    infotokens[_infoinds["fi"]])
            #replace "." with none- missing values. Will turn into null in the db
            valueExtra = [v.strip() if v!="." else None for v in valueExtra]
            value.extend(valueExtra)
#             logger.info(value)  
#             logger.info(len(value))
            values.append(value)
        return(values)
            
    @classmethod
    def getDbWriter(cls,con):
        #building a diplotype writer        
        tablename = "founder_meta"
        cls.buildTable(con, tablename)
        writer = TableWriter(con, cls._outputfields, tablename, INSERT_BATCH_SIZE)
        return(writer)

    @staticmethod
    def buildTable(con, tablename):
        #TODO maxchar
        maxchar = 30
#         "GQ: DP: MQ0F: GP:       AN: MQ: DV: DP4:     SP: SGB:       PV4: FI"
#         15:  4:  0.25: 80,15,0:  2:  25: 4:  0,0,4,0: 0:  -0.556411: .:   0

        con.query("DROP TABLE IF EXISTS " +tablename +";")
        con.query("CREATE TABLE IF NOT EXISTS founder_meta"+
                  "(variant_id integer unsigned NOT NULL, "+ 
                  "strain_name varchar("+str(maxchar)+") NOT NULL, " +
                '''gq smallint unsigned,
                    dp smallint unsigned,
                    mq0f float unsigned,
                    an tinyint unsigned,
                    mq smallint unsigned,
                    dv smallint unsigned,
                    dp_ref_fwd smallint unsigned,
                    dp_ref_rev smallint unsigned,
                    dp_alt_fwd smallint unsigned,
                    dp_alt_rev smallint unsigned,
                    sp smallint unsigned,
                    sgb float, ''' + 
                    #'''PV4 float unsigned  not null,
                    ''' fi tinyint unsigned, 
                    FOREIGN KEY(variant_id) REFERENCES variant(variant_id)
                    ) ENGINE = ''' + engine + ''';''')
    
class InfoParsing:
    _outputFieldToInputField = OrderedDict([("allele",         "Allele"),
                                           ("gene_name",       "Gene"),
                                           ("transcript_name", "Feature"),
                                           ("consequence_name","Consequence")])
    _outputfields = ["variant_id"]
    _outputfields.extend(_outputFieldToInputField.keys())
    
    def __init__(self, headerTokens):
        self.infoCol = headerTokens.index("INFO")
        #This comes out of the comments of the vcf file. Ideally, we would parse it out of the comments rather than hardcoding it here...
        secondaryheader="Allele|Gene|Feature|Feature_type|Consequence|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND"
        secondaryTokens = secondaryheader.split("|")
        #building this map just to use the convenient method for finding header column indexes
        self.colinds = getInputColIndsForOutput(secondaryTokens, InfoParsing._outputFieldToInputField).values()

        self.refCol = headerTokens.index("REF")
        self.altCol = headerTokens.index("ALT")

    
    def getValues(self, vcfTokens, variant_id):
        ref = vcfTokens[self.refCol]
        #the other allele indexes depend on the position within the alt field
        alts = vcfTokens[self.altCol]
        alts = alts.split(",")
        alleles = [ref]
        alleles.extend(alts)
        lengths = [len(allele) for allele in alleles]
        
        minLength  = min(lengths)
        maxLength  = max(lengths)
        indel = minLength!=maxLength
        prefix = ""
        if(indel):
            #vcf reprensents insertion and deletions in the info field with a form of compression, 
            # in which the shortest allele, corresponding to the longest deletion, or no insertion, 
            # is ALWAYS represented in the alleles table as a single base pair, and all other alleles begin with this shortest allele.
            # but in the info field (this is the compression part), the shortest allele is represented as '-',
            # and the other alleles are only implicitly prefixed by the shortest allele. The parsing here will try to remove any implicit 
            # stuff in how our data is stored, and explictly append the prefix to all longer alleles 
            prefix = alleles[lengths.index(minLength)]
        
        info = vcfTokens[self.infoCol]
        info = info.split("CSQ=")[1]
        records = info.split(",")
        values = []
        for csq in records: 
            csq = csq.split("|")  
            value = [variant_id]
            parsedValues = [csq[ind] for ind in self.colinds]
            if parsedValues[0]=="-":
                parsedValues[0]= prefix #replacing the implicit "-" which is just the shortest allele
            else:
                parsedValues[0] = prefix + parsedValues[0] #appending the indel prefix (if its not an indel, this is appending nothing)  
            value.extend(parsedValues)
            values.append(value)
        return(values)
          
    @classmethod      
    def getDbWriter(cls, con):
        tablename = "info"
        cls.buildTable(con, tablename)
        writer = TableWriter(con,  cls._outputfields, tablename, INSERT_BATCH_SIZE)

        return(writer)
    @staticmethod
    def buildTable(con, tablename): 
        con.query("DROP TABLE IF EXISTS " +tablename +";")
        con.query('''CREATE TABLE IF NOT EXISTS ''' + tablename + '''  
                           (variant_id integer unsigned, allele varchar(1000), gene_name varchar(30), transcript_name varchar(60), consequence_name varchar(1000), FOREIGN KEY(variant_id) REFERENCES variant(variant_id)) ENGINE = ''' + engine + ''';''')    
   

class AlleleParsing:
    _outputfields = ["variant_id","allele","allele_index"]
    
    def __init__(self, headerTokens):
        self.refCol = headerTokens.index("REF")
        self.altCol = headerTokens.index("ALT")
    
    def getValues(self, vcfTokens, variant_id):
        
        values = []
        #the reference allele has allele index 0
        value = [variant_id]
        ref = vcfTokens[self.refCol]
        value.append(ref)
        value.append(0)
        values.append(value)
        
        #the other allele indexes depend on the position within the alt field
        alts = vcfTokens[self.altCol]
        alts = alts.split(",")
        for ind in range(0,len(alts)):
            value = [variant_id];
            value.append(alts[ind])
            value.append(ind+1)
            values.append(value)
    
    
        return(values)
          
    @classmethod      
    def getDbWriter(cls, con):
        tablename = "allele"
        cls.buildTable(con, tablename)
        writer = TableWriter(con,  cls._outputfields, tablename, INSERT_BATCH_SIZE)
        return(writer)

    @staticmethod
    def buildTable(con, tablename): 
        con.query("DROP TABLE IF EXISTS " +tablename +";")
        con.query('''CREATE TABLE IF NOT EXISTS ''' + tablename + '''  
                           (variant_id integer unsigned not NULL, 
                           allele varchar(1000) not NULL,  
                           allele_index tinyint unsigned not NULL, 
                           index (allele_index, variant_id),
                           index (variant_id, allele_index),
                           FOREIGN KEY(variant_id) REFERENCES variant(variant_id)) ENGINE = ''' + engine + ''';''')    
      
class TableWriter:    
    def __init__(self, con, fieldnames, tablename, INSERT_BATCH_SIZE):
        self.tablename = tablename
        self.INSERT_BATCH_SIZE = INSERT_BATCH_SIZE
        self.insertQueryRoot = TableWriter._genInsertQueryRoot(fieldnames, tablename)
        self.counter = 0
        self.cur = con.cursor()
        self.con = con
        self.fieldValueBatch = []

    @staticmethod
    def _genInsertQueryRoot(fieldnames, tablename):
        insertVariantsQueryRoot = "INSERT INTO " + tablename + " (" + ", ".join(fieldnames) + ")" + \
            " values " + "(" + ",".join(["%s" for x in range(len(fieldnames))]) + ") "
        return insertVariantsQueryRoot;

    def write(self, fieldvalues):
        
        self.counter = self.counter + 1#len(fieldvalues)
        self.fieldValueBatch.extend(fieldvalues)

        if (self.counter > self.INSERT_BATCH_SIZE):
            self.cur.executemany(self.insertQueryRoot, self.fieldValueBatch)
            self.fieldValueBatch = []
            self.counter = 0
            self.con.commit()
            print("writing batch")
            logger.info("wrote")
    
    def close(self):
        logger.info("closing" + self.tablename)
#         logger.info(self.insertQueryRoot)
#         logger.info(self.fieldValueBatch)
        self.cur.executemany(self.insertQueryRoot, self.fieldValueBatch)
        
        self.fieldValueBatch = []
        self.counter = 0
        self.con.commit()
#         self.con.commit()
        self.cur.close()
    
if __name__=="__main__":    
    logger.info("beginning code")
    parser = argparse.ArgumentParser()
    parser.add_argument('vcfFolder',     help='comma separated vcf files to be broken down')
    parser.add_argument('founderfile',   help='a file with one row per founder name of interest')
    parser.add_argument('db',            help='the database to which we will write the parsed tables')
    parser.add_argument('host',          help="host")
    parser.add_argument('user',          help="user")
    parser.add_argument('passwd',        help="passwd")
    parser.add_argument('-i', '--initial_variant_id', help = '''the starting variant_id, useful if we want to generate unique ids across databases''')
    parser.add_argument('-t', '--allowed_transcripts_file', help ='''an optional file in which each row is the ensembl name of a transcript that will be parsed. 
                                                                     If this argument is included, and a vcf transcript is not in the file, the transcript will be 
                                                                      ignored, and left out of the parsed db.''', action="store", dest="allowed_transcripts_file")
    parser.add_argument('-l', '--limit', help ='''an optional variant count limit; no more than this many variants will be parsed.''', action="store", dest = "limit") 

    parser.add_argument('-e', '--engine', help='''Engine used for all built tables''', action = "store", dest = "engine")
    
    
    args = parser.parse_args()

    logger.info(args)
    ##This is global, TODO, pass it all the way through
    engine = args.engine
    if engine is None:
        engine = "MYISAM" 

    logger.info(args)
    
    limit = args.limit
    if not limit is None:
        limit = int(limit)
    parsevcf(args.vcfFolder, args.founderfile, args.db, args.host, args.user, args.passwd, args.initial_variant_id, args.allowed_transcripts_file,  limit)
    logger.info("done parsing in python code")
