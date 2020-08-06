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
	#Step_6_1: Variant allele frequency caculating
	messageStr = "Step_6_1: Variant allele frequency caculating..."
	ClockAndReport_func(messageStr)
	VafCalculate_func(argElm)
	
	#Step_6_2: RDDpred results report generating
	messageStr = "Step_6_2: RDDpred results report generating..."
	ClockAndReport_func(messageStr)
	ResultReprotGenerate_func(argElm)

	##End main

def ResultReprotGenerate_func(argElm):
	outPrefix = argElm.Out_Prefix_str	
	bamList_path = argElm.RNA_BamList_path
	bamPath_list = open(bamList_path, 'r').read().split()
	vafReport_path = "%s.VafList.txt"%outPrefix
	vafReport_file = open(vafReport_path, 'r')

	sampleList = map(lambda x:"Sample.%s"%(x+1), range(len(bamPath_list)))
	rddReport_path = "%s.RDDpred.results_report.txt"%outPrefix
	rddReport_file = open(rddReport_path, 'w')
	headLine = makeLine(["Chromosome", "One_based_Position", "Reference_Base", "Variant_Base", "Predicted_As", "Prediction_Likelihood_Ratio", "Sample_Frequency", "VAF_Average"] + sampleList, tt)
	rddReport_file.write(headLine + nn)
	while True:
		vafReport_line = vafReport_file.readline()
		if not vafReport_line:
			break
			##End if

		lineElm_list = vafReport_line.split()
		posTag, predLabel, predScore = lineElm_list[:3]
		vafList = map(float, lineElm_list[3:])

		chrId, varPos, editType = posTag.split(".")
		refNuc, varNuc = editType.split(":")

		if predLabel == "Positive":
			predResult = "True_Editing"
		else:
			predResult = "Artefact"
			##End if-else

		sampleCnt = len(filter(lambda x:x>0, vafList))
		vafMean = np.mean(vafList)
		reportLine = makeLine([chrId, varPos, refNuc, varNuc, predResult, predScore, sampleCnt, vafMean] + vafList, tt)  
		rddReport_file.write(reportLine + nn)
		##End while
	##End ResultReprotGenerate_func
	

def VafCalculate_func(argElm):
	outPrefix = argElm.Out_Prefix_str	
	bamList_path = argElm.RNA_BamList_path
	bamPath_list = open(bamList_path, 'r').read().split()
	sampleElm_list = []

	processNum = argElm.Process_Num_int
	VafDir_path = "%s.VafDir"%outPrefix
	makePath(VafDir_path)

	samTools_path = argElm.Tools_Dir_path + ss + "samtools"
	bcfTools_path = argElm.Tools_Dir_path + ss + "bcftools"

	#[1]
	TemplateSiteExtract_sub(argElm)

	#[2]
	VafCalculate_sub(argElm)

	#[3]
	VafSummarize_sub(argElm)

	#[4]
	VafMerge_sub(argElm)
	
	##End VafCalculate_func

def VafMerge_sub(argElm):
	outPrefix = argElm.Out_Prefix_str	
	VafDir_path = "%s.VafDir"%outPrefix
	processNum = int(argElm.Process_Num_int)

	predSiteList_path = VafDir_path + ss + "Predicted.SiteList.txt"
	templateList_path = VafDir_path + ss + "Template.SiteList.txt"

	tempList_path = VafDir_path + ss + "Merged.Template.VafList.txt"
	buffList_path = VafDir_path + ss + "Merged.Buffer.VafList.txt"
	bamList_path = argElm.RNA_BamList_path
	bamPath_list = open(bamList_path, 'r').read().split()

	a, b = templateList_path, tempList_path
	mvCmd = "mv %s %s"%(a, b)
	commands.getoutput(mvCmd)
	for sampleIdx in range(len(bamPath_list)):
		sampleNum = sampleIdx + 1
		sampleVafSite_path = VafDir_path + ss + "Sample.%s.VafList.txt"%sampleNum
		sampleVafList_path = VafDir_path + ss + "Sample.%s.NakedList.txt"%sampleNum

		a, b = sampleVafSite_path, sampleVafList_path
		nakeCmd = "awk '{print $2}' %s > %s"%(a, b)
		commands.getoutput(nakeCmd)
	
		a, b, c = tempList_path, sampleVafList_path, buffList_path
		pasteCmd = "paste %s %s > %s"%(a, b, c)
		commands.getoutput(pasteCmd)

		a, b = buffList_path, tempList_path	
		mvCmd = "mv %s %s"%(a, b)
		commands.getoutput(mvCmd)
		##End for

	predOut_path = "%s.Prediction.ResultList.txt"%outPrefix
	vafReport_path = "%s.VafList.txt"%outPrefix
	a, b, c = predOut_path, tempList_path, vafReport_path
	joinCmd = "join %s %s | tr \" \" \"\t\" > %s"%(a, b, c)
	commands.getoutput(joinCmd)

	rmCmd = "rm -r %s"%VafDir_path
	commands.getoutput(rmCmd)
	##End VafMerge_sub

def VafSummarize_sub(argElm):
	outPrefix = argElm.Out_Prefix_str	
	samTools_path = argElm.Tools_Dir_path + ss + "samtools"
	bcfTools_path = argElm.Tools_Dir_path + ss + "bcftools"
	VafDir_path = "%s.VafDir"%outPrefix
	processNum = int(argElm.Process_Num_int)

	bamList_path = argElm.RNA_BamList_path
	bamPath_list = open(bamList_path, 'r').read().split()
	predSiteList_path = VafDir_path + ss + "Predicted.SiteList.txt"
	templateList_path = VafDir_path + ss + "Template.SiteList.txt"
	a, b = predSiteList_path, templateList_path
	tempCmd = "cut -f3 %s | sort -k1,1 |uniq > %s"%(a, b)
	commands.getoutput(tempCmd)	

	procList = []
	#semElm = mp.Semaphore(processNum)
	semElm = mp.Semaphore(3)
	for sampleIdx in range(len(bamPath_list)):
		sampleNum = sampleIdx + 1
		sampleVaf_path = VafDir_path + ss + "Sample.%s.txt"%sampleNum
		sampleSite_path = VafDir_path + ss + "Sample.%s.SiteList.txt"%sampleNum
		siteCmd = "awk '{print $1}' %s > %s"%(sampleVaf_path, sampleSite_path)
		commands.getoutput(siteCmd)
		
		diffVaf_path = VafDir_path + ss + "Sample.%s.Diff.txt"%sampleNum
		joinVaf_path = VafDir_path + ss + "Sample.%s.Join.txt"%sampleNum
		sampleVafSite_path = VafDir_path + ss + "Sample.%s.VafList.txt"%sampleNum

		a, b, c = templateList_path, sampleVaf_path, joinVaf_path
		d, e, f = templateList_path, sampleSite_path, diffVaf_path
		g, h, i = joinVaf_path, diffVaf_path, sampleVafSite_path
		vafCmd = "join %s %s | awk '{split($2,a,\",\"); sum=a[1]+a[2]+a[3]+a[4]; print $1, (a[3]+a[4])/sum}'> %s; diff %s %s | grep \"<\" | awk '{print $2, 0.0}' > %s; sort -k1,1 %s %s | uniq > %s"%(a, b, c, d, e, f, g, h, i)

		objectFunc = lambda x,y: (commands.getoutput(x), y.release())
		procElm = mp.Process(target=objectFunc, args=[vafCmd, semElm])

		semElm.acquire()
		procElm.start()
		procList.append(procElm)
		##End for

	for procElm in procList:
		procElm.join()
		##End for	
	##End VafSummarize_sub

def VafCalculate_sub(argElm):
	outPrefix = argElm.Out_Prefix_str	
	samTools_path = argElm.Tools_Dir_path + ss + "samtools"
	bcfTools_path = argElm.Tools_Dir_path + ss + "bcftools"
	VafDir_path = "%s.VafDir"%outPrefix
	bamList_path = argElm.RNA_BamList_path
	predSiteList_path = VafDir_path + ss + "Predicted.SiteList.txt"
	processNum = argElm.Process_Num_int
	RefGenome_path = argElm.Ref_Genome_path

	shellScript_path = VafDir_path + ss + "run.VafCalculate.sh"
	shellScript_file = open(shellScript_path, 'w')
	bamPath_list = open(bamList_path, 'r').read().split()
	sampleNum = 1
	for bamPath in bamPath_list:
		outVcf_path = VafDir_path + ss + "Sample.%s.txt"%sampleNum
		a, b, c, d, e, f = samTools_path, predSiteList_path, RefGenome_path, bamPath, bcfTools_path, outVcf_path
		callCmd = "%s mpileup -l %s -AOs -uf %s %s | %s call -vc -p 2 -V indels | grep -v \\\"#\\\" | cut -f1,2,4,5,8 | awk \\\'{split($4,a,\\\",\\\");print $1\\\".\\\"$2\\\".\\\"$3\\\":\\\"a[1], $5}\\\' | awk \\\'{split($2,a,\\\";\\\");for(i=1;i<=length(a);i++) {split(a[i],b,\\\"=\\\");if(b[1]==\\\"DP4\\\") {c=b[2]}}; print $1, c}\\\' | sort -k1,1 |uniq  > %s"%(a, b, c, d, e, f)
		shellScript_file.write(callCmd + nn)
		sampleNum += 1
		##End for
	shellScript_file.close()
	
	a, b = shellScript_path, processNum
	shellCmd = "cat %s | xargs -I token -P %s sh -c \"token\""%(a, b)
	commands.getoutput(shellCmd)
	##End VafCalculate_sub

def TemplateSiteExtract_sub(argElm):
	outPrefix = argElm.Out_Prefix_str	
	bamList_path = argElm.RNA_BamList_path
	bamPath_list = open(bamList_path, 'r').read().split()
	sampleElm_list = []

	processNum = argElm.Process_Num_int
	VafDir_path = "%s.VafDir"%outPrefix
	predOut_path = "%s.Prediction.ResultList.txt"%outPrefix

	predSiteList_path = VafDir_path + ss + "Predicted.SiteList.txt"
	predSiteList_file = open(predSiteList_path, 'w')
	predOut_file = open(predOut_path, 'r')	
	while True:
		predOut_line = predOut_file.readline()
		if not predOut_line:
			break
			##End if

		posTag, predResult, predScore = predOut_line.split()
		chrId, varPos, editType = posTag.split(".")
		siteLine = makeLine([chrId, varPos, posTag, predScore], tt)
		predSiteList_file.write(siteLine + nn)
		##End while

	predSiteList_file.close()
	##End TemplateSiteExtract_sub


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
