##Parses a VCF file containing founders down into smaller pieces, stores the smaller pieces in a
##Database, which is fully normalized an indexed later on.

DEBUG = True

import os
# if(DEBUG):
#     os.chdir("../")
import gzip; 
import time;

import argparse;

from collections import OrderedDict
from collections import defaultdict
import logging.config
import yaml
import ipdb


logging.config.dictConfig(yaml.load(open("../../nd2016/config/logging.config",'r')))
logger=logging.getLogger("mylog")
logger.debug("logging!")
    
INSERT_BATCH_SIZE = 2500

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


def parsevcf(vcfFiles, foundersFile, dbname, firstVariantId, allowedTranscriptsFile, limit):
    logger.debug("!!!!!!!!!!!!!!!!!!!!!")
    logger.debug("about to begin parsing for db " + dbname)
    tmStart = time.time()

    allowedTranscripts = []
    if not allowedTranscriptsFile is None:
        fh = open(allowedTranscriptsFile)
        for transcript in fh:
            allowedTranscripts.append(transcript.strip())

    
    if limit is None:
        limit = float("inf")

    logger.debug("limit is" + str(limit))
    flag = True    
    if flag:
        
        import glob;
        vcfFiles = os.path.join(vcfFiles,"*vcf.gz")
        logger.debug("reading from " + vcfFiles)
        vcfFiles = glob.glob(vcfFiles)
        logger.debug(vcfFiles)
        headerTokensGlobal = getFirstHeader(vcfFiles)
        logger.debug("got header tokens: ") 
        
        diploParser = DiploParsing(headerTokensGlobal, foundersFile)
        parsers = {}
        parsers["Variant"] = VariantParsing(headerTokensGlobal)
        parsers["Allele"]  = AlleleParsing(headerTokensGlobal)
        parsers["Info"]    = InfoParsing(headerTokensGlobal)
        parsers["Diplo"]   = diploParser

        variant_id = 1
        if(not firstVariantId is None):
            variant_id= int(firstVariantId)
            
        logger.debug("about to iterate over vcffiles")
        for vcfFile in vcfFiles:
            logger.debug("parsing " + vcfFile) 
            if variant_id>limit:
                print("past limit")
##                break
            logger.debug(vcfFile)
            vcfFile = gzip.open(vcfFile, 'r')
            headerTokens = navigateToheaderTokens(vcfFile)
            if ",".join(headerTokens)!=",".join(headerTokensGlobal):
                raise("inconsistent headers between vcf files")

            
            # outputFieldToInputField = OrderedDict([("chrom" , "#CHROM"),
            #                                         ("pos"   , "POS")]) 

            # for founder in founders:
            #     outputFieldToInputField[founder] = founder


            
            # colind = getInputColIndsForOutput(headerTokens, _outputFieldToInputField).values()

            
            for line in vcfFile:
                print("parsealine")

                if variant_id>limit:
                    print("past limit")
##                    break
                ##logger.debug(line)
                vcfTokens = line.split("\t")
                if diploParser.isValidVariant(vcfTokens):
##                    diplo.value = diploParser.getValues(vcfTokens, variant_id)

                    values = {}
                    for parserName, parser in parsers.iteritems():
                        values[parserName] = parser.getValues(vcfTokens, variant_id)
                        
#                    ipdb.set_trace()
                        
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
            

  f(gp=="."):
           return probs
alues(self, vcftokens, variant_id):
        vales = []
        forunder, colind in self.founderinds.iteritems():
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
        values = defaultdict(list)
        transcripts = set()
        for csq in records: 
            csq = csq.split("|")  
            parsedValues = [csq[ind] for ind in self.colinds]
            if parsedValues[0]=="-":
                parsedValues[0]= prefix #replacing the implicit "-" which is just the shortest allele
            else:
                parsedValues[0] = prefix + parsedValues[0] #appending the indel prefix (if its not an indel, this is appending nothing)

            ##For every unique transcript, add a reference consequence for the reference allele.
            if (not (parsedValues[2] in transcripts)):
                referenceValues = list(parsedValues)
                referenceValues[0] = ref
                referenceValues[3] = "reference"
                values[parsedValues[0]].append(referenceValues)
                transcripts.add(parsedValues[2])
                
            values[parsedValues[0]].append(parsedValues[1:4])
            

##        transcripts = transcripts.unique()
        if len(transcripts)>1 and len(values.keys())>1:
            ipdb.set_trace()
        
        return(values)

class AlleleParsing:
    _outputfields = ["variant_id","allele","allele_index"]
    
    def __init__(self, headerTokens):
        self.refCol = headerTokens.index("REF")
        self.altCol = headerTokens.index("ALT")
    
    def getValues(self, vcfTokens, variant_id):
        
        values = {}
        #the reference allele has allele index 0
        ref = vcfTokens[self.refCol]
        values[0] = ref
        
        #the other allele indexes depend on the position within the alt field
        alts = vcfTokens[self.altCol]
        alts = alts.split(",")
        for ind in range(0,len(alts)):
            values[ind+1] = alts[ind]
        return(values)
    
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

if DEBUG:
    
    engine = "myisam"
    parsevcf("../../nd2016/data/vcf/rel1410/", "../../nd2016/output/tmp/afounderfile1a1d2f4c7b5e", "imprinted_1410", None,None,1000)
    ipdb.set_trace()
    x=1/0
    
if __name__=="__main__":    
    logger.info("beginning code")
    parser = argparse.ArgumentParser()
    parser.add_argument('vcfFolder',     help='comma separated vcf files to be broken down')
    parser.add_argument('founderfile',   help='a file with one row per founder name of interest')
    parser.add_argument('db',            help='the database to which we will write the parsed tables')
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
    parsevcf(args.vcfFolder, args.founderfile, args.db, args.initial_variant_id, args.allowed_transcripts_file,  limit)
    logger.info("done parsing in python code")
