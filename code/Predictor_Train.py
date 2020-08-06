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
	#Step_4_1: Training predictor
	messageStr = "Step_4_1: Training predictor..."
	ClockAndReport_func(messageStr)
	PredictorTrain_func(argElm)

	#Step_4_2: Prediction target preparation
	messageStr = "Step_4_2: Prediction target preparation..."
	ClockAndReport_func(messageStr)
	PredictionTargetPrepare_func(argElm)

	##End main

def PredictionTargetPrepare_func(argElm):
	outPrefix = argElm.Out_Prefix_str
	modelDir = "%s.ModelDir"%outPrefix	
	predDir = "%s.PredDir"%outPrefix
	makePath(predDir)

	rawList_path = "%s.RDD.RawList.txt"%outPrefix
	targetList_path = predDir + ss + "Prediction.TargetList.txt"

	targetSize_int = PredTargetExtact_sub(rawList_path, targetList_path)
	##End PredictionTargetPrepare_func

FeatureTag_list = ["ReadDepth", "VAF", "CallQual", "FQ", "SGB", "MQ", "MQB", "MQ0F", "BQB", "VDB", "RPB", "PV1", "PV2", "PV3", "PV4"]
def PredTargetExtact_sub(rawList_path, targetList_path):
	rawList_file = open(rawList_path, 'r')
	targetList_file = open(targetList_path, 'w')
	lineCounter = 0

	while True:
		rawLine = rawList_file.readline()
		if not rawLine:
			break
			##End if

		try:
			posTag, samQual, featLine = rawLine.split()
		except:
			continue
			##End try-except	

		featDic = {}
		for featElm in featLine.split(";"):
			try:
				featTag, featVal = featElm.split("=")
			except:
				continue
				##End try-except

			featDic[featTag] = featVal
			##End for

		try:
			featDic["ReadDepth"] = featDic["DP"]
			featDic["CallQual"] = float(samQual)
			a, b, c, d = map(float, featDic["DP4"].split(cc))
			vafVal = (c+d)/(a+b+c+d)
			featDic["VAF"] = vafVal
			a, b, c, d = map(float, featDic["PV4"].split(cc))
			for i in range(1, 5):
				featDic["PV%s"%i] = [a, b, c, d][i-1]
				##End for
			featVal_list = map(lambda x:round_figures(featDic[x]), FeatureTag_list) 
		except:
			continue
			##End try-except
	
		featVal_line = posTag + tt + makeLine(featVal_list, cc)
		targetList_file.write(featVal_line + nn)
		lineCounter += 1
		##End while

	targetList_file.close()
	return lineCounter
	##End PredTargetExtact_sub

def round_figures(x): 
	y = float("%.5f"%float(x))
	return y		
	##End round_figures


def PredictorTrain_func(argElm):
	outPrefix = argElm.Out_Prefix_str
	modelDir = "%s.ModelDir"%outPrefix	
	trainList_path = modelDir + ss + "Train.RddList.csv"
		
	modelObj_path = modelDir + ss + "Predictor.model"
	trainLog_path = "%s.Step4.Predictor.Training.Log"%outPrefix
	wekaJar_path = argElm.Tools_Dir_path + ss + "weka.jar"
	memLimit = argElm.Memory_UseLimit_str.upper()
	if memLimit[-1] != "G":
		memLimit += "G"
		##End if

	a, b, c, = memLimit, wekaJar_path, trainList_path
	d, e = modelObj_path, trainLog_path 
	trainCmd = "java -Xmx%s -cp %s weka.classifiers.trees.RandomForest -i -k -t %s -d %s > %s"%(a, b, c, d, e)
	commands.getoutput(trainCmd)

	attrEvalLog_path = "%s.Step4.Attribute.Evaluation.Log"%outPrefix
	a, b, c, d = memLimit, wekaJar_path, trainList_path, attrEvalLog_path
	evalCmd = "java -Xmx%s -cp %s weka.attributeSelection.InfoGainAttributeEval -i %s > %s"%(a, b, c, d)
	commands.getoutput(evalCmd)
	##End PredictorTrain_func

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
	##End main


