#!/usr/bin/python
#!/usr/local/bin/python

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

def main():
	argElm = OptionParsing_func()
	ArgumentCheck_func(argElm)

	import Prerequisite_Check as step1
	import Raw_RDD_Identify as step2
	import Training_Data_Prepare as step3
	import Predictor_Train as step4
	import RDD_Predict as step5
	import RDDpred_Report as step6

	#Step_1: Prerequisite_Check.py
	messageStr = "Step_1: Prerequisite_Check.py..."
	ClockAndReport_func(messageStr)
	step1.main(argElm)

	#Step_2: Raw_RDD_Identify.py
	messageStr = "Step_2: Raw_RDD_Identify.py..."
	ClockAndReport_func(messageStr)
	step2.main(argElm)

	#Step_3: Training_Data_Prepare.py
	messageStr = "Step_3: Training_Data_Prepare.py..."
	ClockAndReport_func(messageStr)
	step3.main(argElm)

	#Step_4: Predictor_Train.py
	messageStr = "Step_4: Predictor_Train.py..."
	ClockAndReport_func(messageStr)
	step4.main(argElm)
	
	#Step_5: RDD_Predict.py
	messageStr = "Step_5: RDD_Predict.py..."
	ClockAndReport_func(messageStr)
	step5.main(argElm)	

	#Step_6: RDDpred_Report.py
	messageStr = "Step_6: RDDpred_Report.py..."
	ClockAndReport_func(messageStr)
	step6.main(argElm)

	#Closing Announcement
	messageStr = "\nRDDpred: All procedures are completed."
	ClockAndReport_func(messageStr)
	##End main


import argparse
def OptionParsing_func():
	parser = argparse.ArgumentParser()      #initialize argparse
	parser.add_argument('-rbl', '--RNA_Bam_List', dest="RNA_BamList_path",help = "")
	parser.add_argument('-gsl', '--Genomic_SNPs_List', dest="Genomic_SnpList_path",help = "", default=None)
	parser.add_argument('-rsf', '--Reference_Sequence_Fasta', dest="Ref_Genome_path",help = "")
	parser.add_argument('-tdp', '--Tools_Dir_Path', dest="Tools_Dir_path",help = "")
	parser.add_argument('-pni', '--Process_Num', dest="Process_Num_int",help = "", default=1)
	parser.add_argument('-ops', '--Out_Prefix', dest="Out_Prefix_str",help = "")
	parser.add_argument('-psl', '--Positive_Sites_List', dest="Pos_SiteList_path",help = "")
	parser.add_argument('-nsl', '--Negative_Sites_List', dest="Neg_SiteList_path",help = "")
	parser.add_argument('-tml', '--Train_Max_Limit', dest="Train_MaxLimit_int",help = "", default=100000)
	parser.add_argument('-mul', '--Memory_Usage_Limit', dest="Memory_UseLimit_str",help = "", default="10G")
	parser.add_argument('-sud', '--Storage_Usage_Degree', dest="Storage_UseDegree_str",help = "", default=10)

	argElm = parser.parse_args()
	return argElm
	##End OptionParsing_func

def ArgumentCheck_func(argElm):
	errMessage = "Usage: python RDDpred.py [arguments]" + nn
	errMessage += nn + "Essential arguments:" + nn
	errMessage += " -rbl, --RNA_Bam_List [FILE]\t\tList of sorted_bam. Make sure of having identical bam-headers." + nn
	errMessage += " -rsf, --Reference_Sequence_Fasta [FILE]\t\tFasta file used for alignments" + nn
	errMessage += " -tdp, --Tools_Dir_Path [DIR]\t\tTools directory path" + nn
	errMessage += " -ops, --Out_Prefix [STR]\t\tPrefix string for intermediate files (ex: [prefix].RDD.RawList.txt)" + nn
	errMessage += " -psl, --Positive_Sites_List [FILE]\t\tPublished true-sites (see README)" + nn
	errMessage += " -nsl, --Negative_Sites_List [FILE]\t\tCalcualted false-sites (see README)" + nn

	errMessage += nn + "Optional arguments:" + nn
	errMessage += " -gsl, --Genomic_SNPs_List [FILE]\t\tList of Genomic SNPs wanted to be excluded" + nn
	errMessage += " -pni, --Process_Num [INT]\t\tProcess numbers to use (default:1)" + nn
	errMessage += " -tml, --Train_Max_Limit [INT]\t\tTraining data size upper-limits (default:100000)" + nn
	errMessage += " -mul, --Memory_Usage_Limit [GB]\t\tMaximum usage of memory (default:10G)" + nn
	errMessage += " -sud, --Storage_Usage_Degree [FLOAT]\t\tEX: -sud 1 means using 1-bamfile size of storage (default:10)" + nn

	if argElm.RNA_BamList_path == None:
		print errMessage
		sys.exit()
		##End if
	if argElm.Ref_Genome_path == None:
		print errMessage
		sys.exit()
		##End if
	if argElm.Out_Prefix_str == None:
		print errMessage
		sys.exit()
		##End if
	if argElm.Tools_Dir_path == None:
		print errMessage
		sys.exit()
		##End if
	if argElm.Pos_SiteList_path == None:
		print errMessage
		sys.exit()
		##End if
	if argElm.Neg_SiteList_path == None:
		print errMessage
		sys.exit()
		##End if

	argLog_path = "%s.Argument.Log"%argElm.Out_Prefix_str
	argLog_file = open(argLog_path, 'w')
	argLog_file.write("RNA_Bam_List\t%s"%argElm.RNA_BamList_path + nn) 
	argLog_file.write("Genomic_SNPs_List\t%s"%argElm.Genomic_SnpList_path + nn) 
	argLog_file.write("Reference_Sequence_Fasta\t%s"%argElm.Ref_Genome_path + nn) 
	argLog_file.write("Tools_Dir_Path\t%s"%argElm.Tools_Dir_path + nn) 
	argLog_file.write("Process_Num\t%s"%argElm.Process_Num_int + nn) 
	argLog_file.write("Out_Prefix\t%s"%argElm.Out_Prefix_str + nn) 
	argLog_file.write("Positive_Sites_List\t%s"%argElm.Pos_SiteList_path + nn) 
	argLog_file.write("Negative_Sites_List\t%s"%argElm.Neg_SiteList_path + nn) 
	argLog_file.write("Train_Max_Limit\t%s"%argElm.Train_MaxLimit_int + nn) 
	argLog_file.write("Memory_Usage_Limit\t%s"%argElm.Memory_UseLimit_str + nn) 
	argLog_file.write("Storage_Usage_Degree\t%s"%argElm.Storage_UseDegree_str + nn) 

	argLog_file.close()
	##End ArgumentCheck_func


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
	try: os.mkdir(dirPath)
	except: pass
	##End makePath

def makeLine(tokenList, sepToken):
	return sepToken.join(map(str, tokenList))	
	##End makeLine

if __name__ ==  "__main__" :
	main()
	sys.exit()
