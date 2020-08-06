import os
import sys
import commands
import numpy as np
import time
import math
import multiprocessing as mp

nn = "\n"
tt = "\t"
ss = "/"
cc = ","


def main(argElm):
	#Step_1: Prerequisite_Check.py
	#Step_2: Raw_RDD_Idnetify.py
	
	#[2]
	messageStr = "Step_2_1: Extracting every detectable RDDs from alignments..."
	ClockAndReport_func(messageStr)	
	TaskSplit_dic = RawRddExtract_func(argElm)

	#[3]
	messageStr = "Step_2_2: Excluding possible genomic SNPs based on provided input file..."
	ClockAndReport_func(messageStr)	
	GenomicSnpExclude_func(argElm, TaskSplit_dic)
	##End main

def GenomicSnpExclude_func(argElm, TaskSplit_dic):
	#Step_2_2_1
	GenomicSnpParsing_func(argElm, TaskSplit_dic)

	#Step_2_2_2
	GenomicSnpScreen_func(argElm, TaskSplit_dic)

	#Step_2_2_3
	RddMerge_func(argElm, TaskSplit_dic)

	##End GenomicSnpExclude_func

def RddMerge_func(argElm, TaskSplit_dic):
	outPrefix = argElm.Out_Prefix_str	
	rddSplit_dir = "%s.SplitRdd.Dir"%outPrefix; makePath(rddSplit_dir)
	rawRddList_pre = "%s.RDD.RawList.pre"%outPrefix; open(rawRddList_pre, 'w').close()
	pipeCmd_path = "%s.PipeLine.txt"%outPrefix; pipeCmd_file = open(pipeCmd_path, 'w')
	for chrId in sorted(TaskSplit_dic, key=lambda x:TaskSplit_dic[x][1], reverse=True):
		chrLen, chrSplit_num = TaskSplit_dic[chrId]
		a, b = rddSplit_dir, chrId
		chrRddList_path = "%s/%s.RddList.txt"%(a, b); open(chrRddList_path, 'w').close()
		for splitIdx in range(chrSplit_num + 1):
			pipeLine_cmd = RddMergeCommand_func(chrId, chrLen, chrSplit_num, splitIdx, argElm)
			pipeCmd_file.write(pipeLine_cmd + nn)
			##End for

		a, b, c = chrRddList_path, rawRddList_pre, chrRddList_path
		mergeCmd = "cat %s >> %s; rm %s"%(a, b, c)
		pipeCmd_file.write(mergeCmd + nn)
		##End for

	rawRddList_path = "%s.RDD.RawList.txt"%outPrefix; open(rawRddList_path, 'w').close()
	a, b, c = rawRddList_pre, rawRddList_path, rawRddList_pre
	sortCmd = "sort -k1,1 %s | uniq | tr \" \" \"\t\" > %s; rm %s"%(a, b, c)
	pipeCmd_file.write(sortCmd + nn)

	rmPath_list = map(lambda x:"%s.Split%s.Dir"%(outPrefix, x), ["Bam", "Genomic", "Rdd", "Vcf", "Pos"])
	rmvCmd = "rm -r %s"%makeLine(rmPath_list, " ")
	pipeCmd_file.write(rmvCmd + nn)
	pipeCmd_file.close()

	PipeLineRun_func(pipeCmd_path, 1)
	##End RddMerge_func

def RddMergeCommand_func(chrId, chrLen, chrSplit_num, splitIdx, argElm):
	outPrefix = argElm.Out_Prefix_str	
	rddSplit_dir = "%s.SplitRdd.Dir"%outPrefix

	windowSize = chrLen/chrSplit_num
	windSt, windEd = splitIdx * windowSize, (splitIdx+1) * windowSize	
	windEd = min(windEd, chrLen)
	if windEd - windSt == 0:
		return "echo rescue %s:%s-%s"%(chrId, windSt, windEd)
		##End if

	a, b = rddSplit_dir, chrId
	chrRddList_path = "%s/%s.RddList.txt"%(a, b)
	a, b, c = rddSplit_dir, chrId, splitIdx
	splitRdd_path = "%s/%s.%s.RddList.txt"%(a, b, c)
	a, b, c = splitRdd_path, chrRddList_path, splitRdd_path
	mergeRdd_cmd = "cat %s >> %s; rm %s"%(a, b, c)
	
	return mergeRdd_cmd
	##End RddMergeCommand_func

def GenomicSnpScreen_func(argElm, TaskSplit_dic):
	outPrefix = argElm.Out_Prefix_str	
	rddSplit_dir = "%s.SplitRdd.Dir"%outPrefix; makePath(rddSplit_dir)
	pipeCmd_path = "%s.PipeLine.txt"%outPrefix; pipeCmd_file = open(pipeCmd_path, 'w')
	for chrId in sorted(TaskSplit_dic, key=lambda x:TaskSplit_dic[x][1], reverse=True):
		chrLen, chrSplit_num = TaskSplit_dic[chrId]
		for splitIdx in range(chrSplit_num + 1):
			pipeLine_cmd = GenomicScreen_func(chrId, chrLen, chrSplit_num, splitIdx, argElm)
			pipeCmd_file.write(pipeLine_cmd + nn)
			##End for
		##End for
	pipeCmd_file.close()

	processNum = argElm.Process_Num_int
	PipeLineRun_func(pipeCmd_path, processNum)
	##End GenomicSnpScreen_func

def GenomicScreen_func(chrId, chrLen, chrSplit_num, splitIdx, argElm):
	outPrefix = argElm.Out_Prefix_str	
	vcfExtract_dir = "%s.SplitVcf.Dir"%outPrefix
	rddSplit_dir = "%s.SplitRdd.Dir"%outPrefix
	genomicDir_path = "%s.SplitGenomic.Dir"%outPrefix

	windowSize = chrLen/chrSplit_num
	windSt, windEd = splitIdx * windowSize, (splitIdx+1) * windowSize	
	windEd = min(windEd, chrLen)
	if windEd - windSt == 0:
		return "echo rescue %s:%s-%s"%(chrId, windSt, windEd)
		##End if

	a, b, c = vcfExtract_dir, chrId, splitIdx
	splitVcf_path = "%s/%s.%s.vcf"%(a, b, c)	
	a, b, c = rddSplit_dir, chrId, splitIdx
	splitRna_path = "%s/%s.%s.RnaList.txt"%(a, b, c)
	a, b, c, d, e, f = chrId, windSt, windEd, splitVcf_path, splitRna_path, splitVcf_path
	vcfSort_cmd = "awk '{split($5,a,\",\");b=toupper($4);c=toupper(a[1]);  if($1==\"%s\" &&  $2>=%s && $2<%s) print $1\".\"$2\".\"b\":\"c, $6, $8}' %s | sort -k1,1 | uniq  > %s; rm %s"%(a, b, c, d, e, f)

	a, b, c = genomicDir_path, chrId, splitIdx
	splitGenomic_path = "%s/%s.%s.GenomicList.txt"%(a, b, c)
	a, b, c = rddSplit_dir, chrId, splitIdx
	screenRna_path = "%s/%s.%s.ScreenedList.txt"%(a, b, c)
	a, b, c, d = splitRna_path, splitGenomic_path, screenRna_path, splitGenomic_path
	vcfScreen_cmd = "join %s %s > %s; rm %s"%(a, b, c, d)
	
	a, b, c = rddSplit_dir, chrId, splitIdx
	splitRdd_path = "%s/%s.%s.RddList.txt"%(a, b, c)
	a, b, c, d, e = splitRna_path, screenRna_path, splitRdd_path, splitRna_path, screenRna_path
	rddScreen_cmd = "diff %s %s | grep \"<\" | cut -d\" \" -f2- > %s; rm %s %s"%(a, b, c, d, e)

	return makeLine([vcfSort_cmd, vcfScreen_cmd, rddScreen_cmd], ";;")
	##End GenomicScreen_func

def GenomicSnpParsing_func(argElm, TaskSplit_dic):
	outPrefix = argElm.Out_Prefix_str	
	processNum = argElm.Process_Num_int
	genomicDir_path = "%s.SplitGenomic.Dir"%outPrefix; makePath(genomicDir_path)
	genomicList_path = argElm.Genomic_SnpList_path
	if genomicList_path == None:
		rescuePath = "%s.Step2.Rescue.GenomicList.txt"%outPrefix
		open(rescuePath, 'w').close()
		argElm.Genomic_SnpList_path = rescuePath
		##End if

	pipeCmd_path = "%s.PipeLine.txt"%outPrefix; pipeCmd_file = open(pipeCmd_path, 'w')
	for chrId in sorted(TaskSplit_dic, key=lambda x:TaskSplit_dic[x][1], reverse=True):
		chrLen, chrSplit_num = TaskSplit_dic[chrId]
		for splitIdx in range(chrSplit_num + 1):
			pipeLine_cmd = GenomicSplit_func(chrId, chrLen, chrSplit_num, splitIdx, argElm)
			pipeCmd_file.write(pipeLine_cmd + nn)
			##End for
		##End for
	pipeCmd_file.close()

	processNum = argElm.Process_Num_int
	PipeLineRun_func(pipeCmd_path, processNum)
	##End GenomicSnpParsing_func

def GenomicSplit_func(chrId, chrLen, chrSplit_num, splitIdx, argElm):
	outPrefix = argElm.Out_Prefix_str	
	genomicList_path = argElm.Genomic_SnpList_path
	genomicDir_path = "%s.SplitGenomic.Dir"%outPrefix

	a, b, c = genomicDir_path, chrId, splitIdx
	splitGenomic_path = "%s/%s.%s.GenomicList.txt"%(a, b, c)

	windowSize = chrLen/chrSplit_num
	windSt, windEd = splitIdx * windowSize, (splitIdx+1) * windowSize	
	windEd = min(windEd, chrLen)
	if windEd - windSt == 0:
		return "echo rescue %s:%s-%s"%(chrId, windSt, windEd)
		##End if

	a, b, c, d, e = chrId, windSt, windEd, genomicList_path, splitGenomic_path
	chrSplit_cmd = "awk '{if($1==\"%s\" &&  $2>=%s && $2<%s) print $1\".\"$2\".\"$3\":\"$4}' %s | sort -k1,1 | uniq  > %s"%(a, b, c, d, e)
	return chrSplit_cmd
	##End GenomicSplit_func

import multiprocessing as mp
def PipeLineRun_func(pipeCmd_path, processNum):
	semElm = mp.Semaphore(int(processNum))
	procElm_list = []
	pipeCmd_file = open(pipeCmd_path, 'r')
	while True:
		pipeCmd_line = pipeCmd_file.readline()
		if not pipeCmd_line:
			break
			##End if

		semElm.acquire()
		argVect = pipeCmd_line, semElm
		procElm = mp.Process(target=RunCmd_sub, args=argVect)
		procElm.start()
		procElm_list.append(procElm)
		##End while

	pipeCmd_file.close()
	for procElm in procElm_list:
		procElm.join()
		##End for
	
	rmCmd = "rm %s"%pipeCmd_path
	commands.getoutput(rmCmd)
	##End PipeLineRun_func

import subprocess 
def RunCmd_sub(pipeCmd_line, semElm):
	for subCmd in pipeCmd_line.split(";;"):
		cmdLine = subCmd.strip()
		fNull = open(os.devnull, 'w')
		p = subprocess.call(cmdLine, shell=True, stdout=fNull, stderr=subprocess.STDOUT)
		##Edn for

	semElm.release()
	##End RunCmd_sub

def RawRddExtract_func(argElm):
	#Step_2_1_1
	TaskSplit_dic = TaskSplitting_func(argElm)

	#Step_2_1_2
	SplitCalling_func(argElm, TaskSplit_dic)	

	return TaskSplit_dic
	##End RawRddExtract_func

def SplitCalling_func(argElm, TaskSplit_dic):
	outPrefix = argElm.Out_Prefix_str	
	bamExtract_dir = "%s.SplitBam.Dir"%outPrefix; makePath(bamExtract_dir)
	vcfExtract_dir = "%s.SplitVcf.Dir"%outPrefix; makePath(vcfExtract_dir)
	pipeCmd_path = "%s.PipeLine.txt"%outPrefix; pipeCmd_file = open(pipeCmd_path, 'w')
	for chrId in sorted(TaskSplit_dic.keys(), key=lambda x:TaskSplit_dic[x][0], reverse=True):
		chrLen, chrSplit_num = TaskSplit_dic[chrId]
		for splitIdx in range(chrSplit_num + 1):
			pipeLine_cmd = CommandCombine_func(chrId, chrLen, chrSplit_num, splitIdx, argElm)
			pipeCmd_file.write(pipeLine_cmd + nn)
			##End for
		##End for
	pipeCmd_file.close()

	processNum = argElm.Process_Num_int
	PipeLineRun_func(pipeCmd_path, processNum)
	##End SplitCalling_func


def CommandCombine_func(chrId, chrLen, chrSplit_num, splitIdx, argElm):
	outPrefix = argElm.Out_Prefix_str	
	toolsDir_path = argElm.Tools_Dir_path
	bamList_path = argElm.RNA_BamList_path
	refGenome_path = argElm.Ref_Genome_path

	bamExtract_dir = "%s.SplitBam.Dir"%outPrefix
	vcfExtract_dir = "%s.SplitVcf.Dir"%outPrefix
	splitVcf_path = "%s/%s.%s.vcf"%(vcfExtract_dir, chrId, splitIdx)

	samTools_path = "%s/samtools"%toolsDir_path
	bamTools_path = "bamtools"
	bcfTools_path = "%s/bcftools"%toolsDir_path

	windowSize = chrLen/chrSplit_num
	windSt, windEd = splitIdx * windowSize, (splitIdx+1) * windowSize	
	windEd = min(windEd, chrLen)
	if windEd - windSt == 0:
		return "echo rescue %s:%s-%s"%(chrId, windSt, windEd)
		##End if

	posLine = makeLine([chrId, windSt+1, windEd], tt)
	posTag_a = "%s:%s..%s"%(chrId, windSt, windEd)
	posTag_b = "%s:%s-%s"%(chrId, windSt, windEd)
	bamHeader_path = argElm.BamHeader_path

	bamPath_list = argElm.BamPath_list
	bamSplitExtract_dir = "%s/%s.%s.Split.Dir"%(bamExtract_dir, chrId, splitIdx); makePath(bamSplitExtract_dir)
	splitBamList_path = "%s/%s.%s.BamList.txt"%(bamSplitExtract_dir, chrId, splitIdx)
	splitBamList_file = open(splitBamList_path, 'w')
	bamExtractCmd_list = []

	mergedBam_path = "%s/%s.%s.Merged.bam"%(bamExtract_dir, chrId, splitIdx)
	mergedSam_path = "%s/%s.%s.Merged.sam"%(bamExtract_dir, chrId, splitIdx)
	for bamIdx in range(len(bamPath_list)):
		bamPath = bamPath_list[bamIdx]
		if bamIdx == 0:
			a, b, c, d = bamTools_path, bamPath, posTag_a, mergedBam_path
			e, f, g = samTools_path, mergedBam_path, mergedSam_path
			bamExtract_cmd = "%s filter -in %s -region %s -out %s; %s view -h %s > %s"%(a, b, c, d, e, f, g)
			bamExtractCmd_list.append(bamExtract_cmd)
			continue
			##End if

		splitBam_path = "%s/%s.%s.%s.bam"%(bamSplitExtract_dir, chrId, splitIdx, bamIdx)
		a, b, c, d = bamTools_path, bamPath, posTag_a, splitBam_path
		e, f, g = samTools_path, splitBam_path, mergedSam_path
		bamExtract_cmd = "%s filter -in %s -region %s -out %s; %s view %s >> %s"%(a, b, c, d, e, f, g)
		cmdLine = bamExtract_cmd

		splitBamList_file.write(splitBam_path + nn)
		bamExtractCmd_list.append(cmdLine)
		bamIdx += 1
		##End for
	splitBamList_file.close()

	bamExtract_cmd = makeLine(bamExtractCmd_list, ";;")
	preRemove_cmd = "rm -r %s"%bamSplitExtract_dir

	a, b, c, d = samTools_path, mergedSam_path, mergedBam_path, mergedSam_path
	bamConvert_cmd = "%s view -bS %s > %s; rm %s"%(a, b, c, d)
	sortedBam_path = "%s/%s.%s.Sorted.bam"%(bamExtract_dir, chrId, splitIdx)
	a, b, c, d = samTools_path, mergedBam_path, sortedBam_path[:-4], mergedBam_path
	bamSort_cmd = "%s sort %s %s; rm %s"%(a, b, c, d)
	a, b, c, d, e, f = samTools_path, refGenome_path, sortedBam_path, bcfTools_path, splitVcf_path, sortedBam_path
	vcfExtract_cmd = "%s mpileup -R -AOs -uf %s %s | %s call -vc -V indels -p 2 | cut -f-8 > %s; rm %s"%(a, b, c, d, e, f)

	return makeLine([bamExtract_cmd, preRemove_cmd, bamConvert_cmd,  bamSort_cmd, vcfExtract_cmd], ";;")
	##End CommandCombine_func


def TaskSplitting_func(argElm):
	outPrefix = argElm.Out_Prefix_str	
	toolsDir_path = argElm.Tools_Dir_path
	bamList_path = argElm.RNA_BamList_path
	seqFasta_path = argElm.Ref_Genome_path

	samTools_path = "%s/samtools"%toolsDir_path
	bamTools_path = "bamtools"
	bcfTools_path = "%s/bcftools"%toolsDir_path

	bamLocal_path = "%s.Step2.BamList.txt"%outPrefix
	a, b = bamList_path, bamLocal_path
	commands.getoutput("cp %s %s"%(a, b))
	argElm.RNA_BamList_path = bamLocal_path
	bamList_file = open(bamList_path, 'r'); argElm.BamPath_list = bamList_file.read().split(); bamList_file.close()

	nullSam_path = "%s.Step2.Null.sam"%outPrefix; open(nullSam_path, 'w').close()
	bamHeader_path = "%s.Step2.BamHeader.txt"%outPrefix

	seqMirror_path = "%s.Step2.SeqMirror.fasta"%outPrefix
	a, b = os.path.abspath(seqFasta_path), seqMirror_path	
	seqMirror_cmd = "ln -s %s %s"%(a, b); commands.getoutput(seqMirror_cmd)
	a, b = samTools_path, seqMirror_path
	faIdx_cmd = "%s faidx %s; rm %s"%(a, b, b); commands.getoutput(faIdx_cmd)

	seqIndex_path = "%s.Step2.SeqMirror.fasta.fai"%outPrefix
	a, b, c, d = samTools_path, seqIndex_path, nullSam_path, bamHeader_path
	headerGen_cmd = "%s view -ht %s %s > %s; rm %s %s"%(a, b, c, d, b, c)
	commands.getoutput(headerGen_cmd)
	argElm.BamHeader_path = bamHeader_path

	TaskSplit_dic = {}
	bamHeader_file = open(bamHeader_path, 'r')
	while True:
		bamHeader_line = bamHeader_file.readline()
		if not bamHeader_line:
			break
			##End if

		try: headTag, chrTag, lenTag = bamHeader_line.split()
		except:	continue
		if headTag != "@SQ":
			continue
			##End if

		chrId, chrLen = map(lambda x:x.split(":")[-1], [chrTag, lenTag])
		TaskSplit_dic[chrId] = int(chrLen)
		##End while
	bamHeader_file.close()

	chrLen_sum = sum(map(lambda x:float(TaskSplit_dic[x]), TaskSplit_dic.keys()))
	processNum = int(argElm.Process_Num_int)
	sudDegree = float(argElm.Storage_UseDegree_str)
	bamList_file = open(bamList_path, 'r')
	bamNum = len(bamList_file.read().split())
	bamList_file.close()
	for chrId in TaskSplit_dic:
		chrLen = TaskSplit_dic[chrId]
		chrSize_frac = chrLen/chrLen_sum
		"""
		#chrSize_frac/chrSplit_num * bamNum * processNum < sudDegree
		chrSize_frac/chrSplit_num * bamNum * processNum < sudDegree/50
		(chrSize_frac * bamNum * processNum)/sudDegree < chrSplit_num
		"""
		chrSplit_num = int((chrSize_frac * bamNum * processNum * 50)/sudDegree)
		if not chrSplit_num:
			chrSplit_num = 1
			##End if
		TaskSplit_dic[chrId] = chrLen, chrSplit_num
		##End for

	return TaskSplit_dic
	##End TaskSplitting_func


import datetime as dt
def ClockAndReport_func(messageStr):
	currentTime = dt.datetime.now()
	print "%s [%s]"%(messageStr, currentTime)
	##End ClockAndReport_func


def GetNumString_func(sampleNum, sampleSize):
	decNum = int(math.log(sampleSize, 10)) + 1
	numString = str(sampleNum)
	while len(numString) < decNum:
		numString = "0" + numString
		##End while

	return numString
	##End GetNumString_func

def makePath(dirPath):
	return commands.getoutput("mkdir %s"%dirPath)
	##End makePath

def makeLine(tokenList, sepToken):
	return sepToken.join(map(str, tokenList))	
	##End makeLine

if __name__ ==  "__main__" :
	main()
	sys.exit()
