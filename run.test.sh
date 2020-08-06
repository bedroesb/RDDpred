mkdir -p OutDir
python RDDpred.py -rbl Test.BamList.txt  -rsf hg19_all.fa -tdp ToolBox.Dir -ops OutDir/Test -psl TestBam.Dir/Test.Positive.txt -nsl TestBam.Dir/Test.Negative.txt -pni 6
