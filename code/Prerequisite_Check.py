import os
import sys
import commands
import numpy as np
import time
import math

import multiprocessing as mp
import itertools

nn = "\n"
tt = "\t"
ss = "/"
cc = ","

def main(argElm):
	#Step_1_1: Checking Samtools pipeline
	messageStr = "Step_1_1: Checking Samtools..."
	ClockAndReport_func(messageStr)
	SamtoolsCheck_func(argElm)

	#Step_1_2: Checking Bcftools pipeline
	messageStr = "Step_1_2: Checking Bcftools..."
	ClockAndReport_func(messageStr)
	BcftoolsCheck_func(argElm)

	#Step_1_3: Checking Bamtools pipeline
	messageStr = "Step_1_3: Checking Bamtools..."
	ClockAndReport_func(messageStr)
	BamtoolsCheck_func(argElm)

	#Step_1_4: Checking WEKA pipeline
	messageStr = "Step_1_4: Checking WEKA..."
	ClockAndReport_func(messageStr)
	WekaCheck_func(argElm)

	##End main


#WebPage_address = "biohealth.snu.ac.kr/software/RDDpred"
WebPage_address = "README"
def WekaCheck_func(argElm):
	wekaJar_path = argElm.Tools_Dir_path + ss + "weka.jar"
	testData_dir = argElm.Tools_Dir_path + ss + "TestData.Dir"
	outPrefix = argElm.Out_Prefix_str	
	testResult_dir = "%s.TestResult.Dir"%outPrefix

	testData_path = testData_dir + ss + "Weka.TestData.csv"
	testAnswer_path = testData_dir + ss + "Weka.TestAnswer.txt"
	testResult_path = testResult_dir + ss + "Weka.TestResult.txt"

	correctAnswer_md5sum = commands.getoutput("md5sum %s"%testAnswer_path).split()[0]
	a, b, c = wekaJar_path, testData_path, testResult_path
	commands.getoutput("java -cp %s weka.classifiers.trees.RandomForest -t %s -p 0 > %s"%(a, b, c))
	testResult_md5sum = commands.getoutput("md5sum %s"%testResult_path).split()[0]
	
	if correctAnswer_md5sum != testResult_md5sum:
		errMessage = " [Error_4] Weka is not working properly... (see %s)"%WebPage_address
		sys.exit(errMessage)
		##End if

	removeCmd = "rm -rf %s"%testResult_dir
	commands.getoutput(removeCmd)
	##End SamBcftoolsCheck_func

def BamtoolsCheck_func(argElm):
	#bamTools_path = argElm.Tools_Dir_path + ss + "bamtools"
	bamTools_path = "bamtools"
	testData_dir = argElm.Tools_Dir_path + ss + "TestData.Dir"
	outPrefix = argElm.Out_Prefix_str	
	testResult_dir = "%s.TestResult.Dir"%outPrefix

	testData_path = testData_dir + ss + "Bamtools.TestData.bam"
	testAnswer_path = testData_dir + ss + "Bamtools.TestAnswer.txt"
	testResult_path = testResult_dir + ss + "Bamtools.TestResult.txt"

	correctAnswer_md5sum = commands.getoutput("md5sum %s"%testAnswer_path).split()[0]
	a, b, c = bamTools_path, testData_path, testResult_path
	commands.getoutput("%s header -in %s > %s"%(a, b, c))
	testResult_md5sum = commands.getoutput("md5sum %s"%testResult_path).split()[0]
	
	if correctAnswer_md5sum != testResult_md5sum:
		errMessage = " [Error_3] Bamtools is not working properly... (see %s)"%WebPage_address
		sys.exit(errMessage)
		##End if
	##End SamBcftoolsCheck_func

def BcftoolsCheck_func(argElm):
	bcfTools_path = argElm.Tools_Dir_path + ss + "bcftools"
	testData_dir = argElm.Tools_Dir_path + ss + "TestData.Dir"
	outPrefix = argElm.Out_Prefix_str	
	testResult_dir = "%s.TestResult.Dir"%outPrefix

	testData_path = testData_dir + ss + "Bcftools.TestData.txt"
	testAnswer_path = testData_dir + ss + "Bcftools.TestAnswer.txt"
	testResult_path = testResult_dir + ss + "Bcftools.TestResult.txt"

	correctAnswer_md5sum = commands.getoutput("md5sum %s"%testAnswer_path).split()[0]
	a, b, c = bcfTools_path, testData_path, testResult_path
	commands.getoutput("%s call -c %s | grep -v \"#\" > %s"%(a, b, c))
	testResult_md5sum = commands.getoutput("md5sum %s"%testResult_path).split()[0]
	
	if correctAnswer_md5sum != testResult_md5sum:
		errMessage = " [Error_2] Bcftools is not working properly... (see %s)"%WebPage_address
		sys.exit(errMessage)
		##End if
	##End SamBcftoolsCheck_func

def SamtoolsCheck_func(argElm):
	samTools_path = argElm.Tools_Dir_path + ss + "samtools"
	testData_dir = argElm.Tools_Dir_path + ss + "TestData.Dir"
	outPrefix = argElm.Out_Prefix_str	
	testResult_dir = "%s.TestResult.Dir"%outPrefix
	makePath(testResult_dir)

	testData_path = testData_dir + ss + "Samtools.TestData.bam"
	testAnswer_path = testData_dir + ss + "Samtools.TestAnswer.txt"
	testResult_path = testResult_dir + ss + "Samtools.TestResult.txt"

	correctAnswer_md5sum = commands.getoutput("md5sum %s"%testAnswer_path).split()[0]
	a, b, c = samTools_path, testData_path, testResult_path
	commands.getoutput("%s mpileup %s > %s"%(a, b, c))
	testResult_md5sum = commands.getoutput("md5sum %s"%testResult_path).split()[0]
	
	if correctAnswer_md5sum != testResult_md5sum:
		errMessage = " [Error_1] Samtools is not working properly... (see %s)"%WebPage_address
		sys.exit(errMessage)
		##End if
	##End SamBcftoolsCheck_func

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
