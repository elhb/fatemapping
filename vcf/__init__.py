import sys
sys.stderr.write('module loaded\n')

class VcfEntry(object):
    """Class representing a vcf entry or genomic variant recorded in a vcf file"""    

    def __init__(self,line,sampleNames):

	import time
	import sys

	self.isIndel = None
	self.vcfLine = line
	line = self.vcfLine.rstrip().split('\t')
	self.perSampleInfo = {sampleName:{} for sampleName in sampleNames}

	for i in range(len(line)):
	    if i == 0:   self.chrom = line[i]
	    elif i == 1: self.pos = line[i]
	    elif i == 2:
		self.id = line[i]
		if line[i] == '.':pass
		    #AnalysisPipe.nonDbSNPVarCount+=1;
		    #self .id = 'novel#'+str(AnalysisPipe.nonDbSNPVarCount)
	    elif i == 3: self.refBase = line[i]
	    elif i == 4: self.altBases = line[i].split(',')
	    elif i == 5: self.varQual = line[i]
	    elif i == 6: self.passFilter = line[i]
	    elif i == 7: self.info = line[i]
	    elif i == 8: self.perSampleFormat = line[i]
	    elif i >= 9:
		sampleInfo = line[i]
		infoTags = self.perSampleFormat.split(':')
		if infoTags[0] == 'GT' and sampleInfo[0:2] != './.' or sampleInfo != './.':
		    sampleInfoList = sampleInfo.split(':') #GT:AD:DP:GQ:PL
		    #if sampleInfoList[0] == './.':continue
		    for i2 in range(len(sampleInfoList)): self.perSampleInfo[sampleNames[i-9]][infoTags[i2]] = sampleInfoList[i2]
		    try:
			if self.perSampleInfo[sampleNames[i-9]]['DP'] == '.': self.perSampleInfo[sampleNames[i-9]]['DP'] = 0
			self.perSampleInfo[sampleNames[i-9]]['DP'] = int(self.perSampleInfo[sampleNames[i-9]]['DP'])
		    except KeyError: self.perSampleInfo[sampleNames[i-9]]['DP'] = 0
		    try:
			if self.perSampleInfo[sampleNames[i-9]]['GQ'] == '.': self.perSampleInfo[sampleNames[i-9]]['GQ'] = 0
			self.perSampleInfo[sampleNames[i-9]]['GQ'] = int(self.perSampleInfo[sampleNames[i-9]]['GQ'])
		    except KeyError: self.perSampleInfo[sampleNames[i-9]]['GQ'] = 0
		    if 'GT' not in self.perSampleInfo[sampleNames[i-9]]: self.perSampleInfo[sampleNames[i-9]]['GT'] = 'Unknown'
		    if 'AD' not in self.perSampleInfo[sampleNames[i-9]]: self.perSampleInfo[sampleNames[i-9]]['AD'] = 'Unknown'
		    try:
			if self.perSampleInfo[sampleNames[i-9]]['DP'] == '.': self.perSampleInfo[sampleNames[i-9]]['DP'] = 0
		    except KeyError:
			msg = '#ERROR_MSG#'+time.strftime("%Y-%m-%d:%H:%M:%S",time.localtime())+'# VcfEntry format error: '+str(infoTags)+' '+str(sampleInfoList)+' 
'+str(self.perSampleInfo[sampleNames[i-9]])+' '+str(i)+' >'+str(sampleInfo)+'<'+'.\n'
			sys.stderr.write(msg)

	for altBase in list(self.altBases):
	    if len(self.refBase) != len(altBase):
		self.type = 'INDEL'
		self.isIndel = True
	    elif not self.isIndel:
		self.type = 'SNV'
		self.isIndel = False

def vcfParser(infilename,subsetSize=None,logfile=False):
	"""Generator taking an infilename ie vcf and yielding VcfEntries"""
	
	import os
	import sys
	import time
	from misc import bufcount
	from misc import Progress
	
	nonDbSNPVarCount = 0
	vcfLineCount = 0
	toQcounter = 0
	
	if logfile: logfile.write('#LOGMSG#'+time.strftime("%Y-%m-%d:%H:%M:%S",time.localtime())+'# Loading vcf to memory, my id is '+str(os.getpid())+'.\n')
	lastChomosomeRead = None
	
	if not subsetSize:
	    totalLineNum = bufcount(infilename)
	    if logfile:
		logfile.write('#LOGMSG#'+time.strftime("%Y-%m-%d:%H:%M:%S",time.localtime())+'# '+str(totalLineNum)+' lines found in vcf file.\n')
	    	progress = Progress(totalLineNum, verb='full', logfile=logfile, unit='variants read from disk' ,mem=False, printint=1)
	else:
	    if logfile: progress = Progress(subsetSize, verb='full', logfile=logfile, unit='variants read from disk' ,mem=False)

	with open(infilename) as infile:
	    for line in infile:
		if not subsetSize and logfile: progress.update()
		if subsetSize and toQcounter >= subsetSize:
		    if logfile: logfile.write('#WARNING#'+time.strftime("%Y-%m-%d:%H:%M:%S",time.localtime())+'# Warning running on a subset of variants debbugging purpose.\n')
		    break

                if line[0] == "#": #THIS LINE IS HEADER
                    if line[1] == 'C': #THIS LINE IS HEADER WITH SAMPLE INFO
                        line = line.rstrip().split('\t')
                        sampleNames = [line[i] for i in range(9,len(line))]
                    continue

		vcfLineCount +=1

		variation = VcfEntry(line,sampleNames)

		if variation.passFilter == 'PASS' and not variation.isIndel and variation.varQual >= 30:
		    yield variation		    
		    if subsetSize and logfile: progress.update()
		    toQcounter+=1
		    if variation.chrom != lastChomosomeRead and logfile: logfile.write('#LOGMSG#'+time.strftime("%Y-%m-%d:%H:%M:%S",time.localtime())+'# started loading variations from 
'+str(variation.chrom)+'.\n')
		    lastChomosomeRead = variation.chrom
	if logfile: logfile.write('#LOGMSG#'+time.strftime("%Y-%m-%d:%H:%M:%S",time.localtime())+'# vcf file loaded.\n')
