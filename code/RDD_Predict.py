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
	#Step_5_1: Predicting RDD Candidates
	messageStr = "Step_5_1: Predicting RDD Candidates..."
	ClockAndReport_func(messageStr)
	PredictionPerfom_func(argElm)

	#Step_5_2: Summarizing prediction results
	messageStr = "Step_5_2: Summarizing prediction results..."
	ClockAndReport_func(messageStr)
	PredictionResultSummarize_func(argElm)

	##End main

def PredictionResultSummarize_func(argElm):
	outPrefix = argElm.Out_Prefix_str
	rawRddList_path = "%s.RDD.RawList.txt"%outPrefix
	predDir = "%s.PredDir"%outPrefix
	targetList_path = predDir + ss + "Prediction.TargetList.txt"
	predOut_path = "%s.Prediction.ResultList.txt"%outPrefix

	summaryLog_path = "%s.Step5.Results.Summary.txt"%outPrefix
	summaryLog_file = open(summaryLog_path, 'w')

	rawRdd_count = commands.getoutput("wc -l %s"%rawRddList_path).split()[0]
	summaryLog_file.write("Raw-RDDs:\t%s"%rawRdd_count + nn)
	
	rddTarget_count = commands.getoutput("wc -l %s"%targetList_path).split()[0]
	summaryLog_file.write("Prediction-targets:(>2 variant reads)\t%s"%rddTarget_count + nn)
	
	posPred_count = commands.getoutput("grep -c Positive %s"%predOut_path).strip()
	negPred_count = commands.getoutput("grep -c Negative %s"%predOut_path).strip()
	summaryLog_file.write("Predicted as true-editing(Positive):\t%s"%posPred_count + nn)
	summaryLog_file.write("Predicted as artefacts(Negative):\t%s"%negPred_count + nn)

	try: accRate = float(posPred_count)/float(rawRdd_count) * 100
	except: accRate = float("nan")
	summaryLog_file.write("Overall acceptance rate:\t%.2f%%"%accRate + nn)
	summaryLog_file.close()

	removeCmd = "rm -rf %s"%predDir
	commands.getoutput(removeCmd)
	##End PredictionResultSummarize_func

def PredictionPerfom_func(argElm):
	outPrefix = argElm.Out_Prefix_str
	modelDir = "%s.ModelDir"%outPrefix	
	predDir = "%s.PredDir"%outPrefix

	targetList_path = predDir + ss + "Prediction.TargetList.txt"
	predOut_path = "%s.Prediction.ResultList.txt"%outPrefix
	lineCouter = TargetPredict_sub(argElm)

	#removeCmd = "rm -rf %s"%modelDir
	#commands.getoutput(removeCmd)
	##End PredictionPerfom_func

def TargetPredict_sub(argElm):
	outPrefix = argElm.Out_Prefix_str
	modelDir = "%s.ModelDir"%outPrefix	
	predDir = "%s.PredDir"%outPrefix

	targetList_path = predDir + ss + "Prediction.TargetList.txt"
	predOut_path = "%s.Prediction.ResultList.txt"%outPrefix
	predOut_file = open(predOut_path, 'w')

	lineCounter = 0
	lineBuffer_pair = [], []
	bufferState = "Initiated"
	targetFile = open(targetList_path, 'r')
	while True:
		targetLine = targetFile.readline()
		if not targetLine:
			break
			##End if

		posTag, featLine = targetLine.split()
		if bufferState == "Initiated":
			csvLine = featLine + cc + "Positive"
			lineBuffer_pair[0].append(posTag)
			lineBuffer_pair[1].append(csvLine+nn)
			bufferState = "Activated"
			##End if
		elif bufferState == "Activated":
			csvLine = featLine + cc + "Negative"
			lineBuffer_pair[0].append(posTag)
			lineBuffer_pair[1].append(csvLine+nn)
			##End if

		if len(lineBuffer_pair[0]) == 10**5:
			LineBufferFlush_sub(argElm, lineBuffer_pair, predOut_file)
			lineBuffer_pair = [], []
			bufferState = "Initiated"
			##End if
		##End while

	if len(lineBuffer_pair[0]) > 0:
		LineBufferFlush_sub(argElm, lineBuffer_pair, predOut_file)
		lineBuffer_pair = [], []
		bufferState = "Initiated"
		##End if

	predOut_file.close()
	##End TargetPredict_sub

FeatureTag_list = ["ReadDepth", "VAF", "CallQual", "FQ", "SGB", "MQ", "MQB", "MQ0F", "BQB", "VDB", "RPB", "PV1", "PV2", "PV3", "PV4"]
def LineBufferFlush_sub(argElm, lineBuffer_pair, predOut_file):
	outPrefix = argElm.Out_Prefix_str
	modelDir = "%s.ModelDir"%outPrefix	
	predDir = "%s.PredDir"%outPrefix

	featTag_line = makeLine(FeatureTag_list + ["TrainLabel"], cc) 
	bufferFile_path = predDir + ss + "TargetLine.Buffer.csv"
	bufferFile = open(bufferFile_path, 'w')
	bufferFile.write(featTag_line + nn)
	bufferFile.writelines(lineBuffer_pair[1])
	bufferFile.close()

	modelObj_path = modelDir + ss + "Predictor.model"
	wekaJar_path = argElm.Tools_Dir_path + ss + "weka.jar"
	predOut_path = "%s.Prediction.ResultList.txt"%outPrefix

	a, b, c, d = argElm.Memory_UseLimit_str, wekaJar_path, bufferFile_path, modelObj_path
	predCmd = "java -Xmx%s -cp %s weka.classifiers.trees.RandomForest -T %s -l %s -p 0"%(a, b, c, d)
	predLine_list = commands.getoutput(predCmd).split(nn)
	for predLine in predLine_list:
		lineElm_list = predLine.split()
		try:
			predIdx, x, predLabel = lineElm_list[:3]
			predIdx = int(predIdx) - 1
			predLabel = predLabel.split(":")[1]
			predScore = lineElm_list[-1]
			posTag = lineBuffer_pair[0][int(predIdx)]
		except:
			continue
			##End try-except

		chrId, varPos, editType = posTag.split(".")
		refNuc, varNuc = editType.split(":")

		#predOut_line = makeLine([chrId, varPos, refNuc, varNuc, predLabel, predScore], tt)
		predOut_line = makeLine([posTag, predLabel, predScore], tt)
		predOut_file.write(predOut_line + nn)
		##End for
	##End LineBufferFlush_sub

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
