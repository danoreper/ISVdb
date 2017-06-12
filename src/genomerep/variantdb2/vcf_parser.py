##Parses a VCF file containing founders down into smaller pieces, stores the smaller pieces in a
##Database, which is fully normalized an indexed later on.

DEBUG = False

import readline
import rlcompleter
readline.parse_and_bind("tab: complete")

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
##import ipdb
import cProfile;

logging.config.dictConfig(yaml.load(open("../../nd2016/config/logging.config",'r')))
logger=logging.getLogger("mylog")
logger.debug("logging!")
    
INSERT_BATCH_SIZE = 2500000

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
    logger.debug("trying to pull out header tokens. current path is: " + os.getcwd())
    logger.debug("opening:" + vcfFiles[0])
    ##grab the first vcf, we assume all vcfs will have the same header tokens. Will check this later on
    fle = gzip.open(vcfFiles[0], 'r')
    logger.debug("opened file")
    headerTokensGlobal = navigateToheaderTokens(fle)
    fle.close()
    return headerTokensGlobal


def mergeAndWrite(values, writer):
    ##merge the parsed values together
    for diplo in values["Diplo"]:
        allele_index_1 = diplo[1]
        allele_index_2 = diplo[2]

        if allele_index_1!=None:
            allele_1 = values["Allele"][allele_index_1]
            consequence_1s = values["Info"][allele_1]
        else:
            allele_1 = "None"
            consequence_1s = values["Info"].values()[0]

        if allele_index_2!=None:
            allele_2 = values["Allele"][allele_index_2]
            consequence_2s = values["Info"][allele_2]
        else:
            allele_2 = "None"
            consequence_2s = values["Info"].values()[0]

        for c1, c2 in zip(consequence_1s, consequence_2s):
            if c1[1]!=c2[1]:
                raise("invalid assumption about transcript ordering")

            conseq1 = str(c1[2])
            conseq2 = str(c2[2])
            
            if allele_index_1 == None:
                conseq1 = "None"
                if diplo[0]=="C57BL6J":
                   conseq1 = "reference"

            if allele_index_2 == None:
                conseq2 = "None"
                if diplo[0]=="C57BL6J":
                   conseq2 = "reference"
                   
            prob = "%.3f" %(diplo[3])
            ##prob = str(diplo[3])
            collatedrow = values["Variant"] + [diplo[0], allele_1, allele_2, prob] + c1[0:2] + [conseq1, conseq2]
            # if diplo[0]=="C57BL6J":
            #     ipdb.set_trace()
            #     x = 4
            #     print(5)
            
            collatedrow = "\t".join(collatedrow)+"\n"
            writer.write(collatedrow)


def parsevcf(vcfFiles, foundersFile, dbname, allowedTranscriptsFile, limit):
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
        # logger.debug("got diplo parser")
        parsers = {}
        parsers["Variant"] = VariantParsing(headerTokensGlobal)
        parsers["Allele"]  = AlleleParsing(headerTokensGlobal)
        parsers["Info"]    = InfoParsing(headerTokensGlobal)
        parsers["Diplo"]   = diploParser

        variant_id = 1
        logger.debug("about to iterate over vcffiles")
        writer = open(dbname, "w", INSERT_BATCH_SIZE)
        writer.write("\t".join(["variant_id", "pos", "strain", "allele_1", "allele_2", "prob", "gene_name", "transcript_name", "consequence_1", "consequence_2\n"]))

        for vcfFile in vcfFiles:
            logger.debug("parsing " + vcfFile) 
            if variant_id>limit:
                print("past limit")
                break
            logger.debug(vcfFile)
            vcfFile = gzip.open(vcfFile, 'r')
            headerTokens = navigateToheaderTokens(vcfFile)
            if ",".join(headerTokens)!=",".join(headerTokensGlobal):
                raise("inconsistent headers between vcf files")


            for line in vcfFile:

                if variant_id>limit:
                    print("past limit")
                    break
                
                ##logger.debug(line)
                vcfTokens = line.split("\t")
                if diploParser.isValidVariant(vcfTokens):
                    # logger.debug("valid line!")
                    # logger.debug(line)
                    # logger.debug(vcfTokens)
                    values = {}
                    for parserName, parser in parsers.iteritems():
                        values[parserName] = parser.getValues(vcfTokens, variant_id)
                    mergeAndWrite(values, writer)

                variant_id = variant_id + 1
                if (variant_id % 10000)==0:
                    print(variant_id)
            vcfFile.close()
        writer.close()
            
        #close writers only after we are done with all vcfs
    logger.debug("attempting to postprocess")
    tmEnd = time.time();
    logger.debug(tmEnd-tmStart)

class VariantParsing:    
    _outputFieldToInputField = OrderedDict([("pos"   , "POS")])  
    _outputfields = ["variant_id"]
    _outputfields.extend(_outputFieldToInputField.keys())
     

    def __init__(self, headerTokens):      
        self.colind = getInputColIndsForOutput(headerTokens, VariantParsing._outputFieldToInputField).values()
        
    def getValues(self, vcftokens, variant_id):
        value = [str(variant_id)];
        value.extend([vcftokens[ind] for ind in self.colind])
        return(value)

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
##        ipdb.set_trace()

        
        for colind in self.founderinds.itervalues():
        ##    diplotypeVal = vcftokens[colind].split(":")
        ##    diplotypeVal = diplotypeVal[0]
            tok = vcftokens[colind]
            if not tok.startswith("0/0") and not tok.startswith("./."):
                return(True)

        return False


    def getProbs(self, gps):
        probs = len(gps) * [0]
        denom = -.1
        for i in range(0, len(gps)):
            gp = gps[i]
            if(gp=="."):
                probs[i] = 0
            else:
                probs[i] = pow(10, int(gp) * denom)
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

            if mleDiplotype=="./.":
                values.append([founder, None, None, 1])
            elif all([x == "." for x in gps]):
                values.append([founder, None, None, 1])

            else:
                probs = self.getProbs(gps)
                probSum = sum(probs)
                denom = 1.0/probSum
#                probs = [prob*denom for prob in probs]

                for i in range(0,len(gps)):
                    prob = probs[i]*denom
                    if prob>=.001:
                        allele_index_1 = self.jalleles_by_F[i]
                        allele_index_2 = self.kalleles_by_F[i]
                        values.append([founder, allele_index_1, allele_index_2, prob])
        values.append(["C57BL6J", 0,0,1])
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
##                if(ref =="Info"):
##                    ipdb.set_trace()
##                ipdb.set_trace()
                values[ref].append(referenceValues[1:4])
                transcripts.add(parsedValues[2])
                
            values[parsedValues[0]].append(parsedValues[1:4])

##        transcripts = transcripts.unique()
##        if len(transcripts)>1 and len(values.keys())>1:
##            ipdb.set_trace()
        
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


if DEBUG:
    parsevcf("../../nd2016/data/isvdb/exon_1410/temp/vcf/chr_1",
             "../../nd2016/data/isvdb/exon_1410/temp/vcf/chr_1/founder.csv",
             "../../nd2016/data/isvdb/exon_1410/temp/founders/chr_1.txt",
             None,
             100000)
##    ipdb.set_trace()
    x=1/0
    
if __name__=="__main__":    
    logger.info("beginning code")
    parser = argparse.ArgumentParser()
    parser.add_argument('vcfFolder',     help='comma separated vcf files to be broken down')
    parser.add_argument('founderfile',   help='a file with one row per founder name of interest')
    parser.add_argument('db',            help='the database to which we will write the parsed tables')
    
    parser.add_argument('-t', '--allowed_transcripts_file',
                        help ='''an optional file in which each row is the ensembl name of a transcript that will be parsed. If this argument is included, and a vcf transcript is not in the file, the transcript will be 
    ignored, and left out of the parsed db.''', action="store", dest="allowed_transcripts_file")
    
    parser.add_argument('-l', '--limit', help ='''an optional variant count limit; no more than this many variants will be parsed.''', action="store", dest = "limit") 
        
    args = parser.parse_args()

    logger.info(args)
    ##This is global, TODO, pass it all the way through
    logger.info(args)
    
    limit = args.limit
    
    if not limit is None:
        limit = int(limit)
    parsevcf(args.vcfFolder, args.founderfile, args.db, args.allowed_transcripts_file,  limit)
    logger.info("done parsing in python code")
