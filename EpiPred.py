#!/usr/bin/env python

# Author:
#    Frederick K Kamanu, Ph.D [frederick dot kamanu at gmail.com]
# Contributors:
#   Geoffrey Siwo, Ph.D
# All Rights Reserved.

#*****************************************************************************************
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/>.
#*****************************************************************************************


import sys
import os
import re
import subprocess
import glob
from optparse import OptionParser


#Main Class
#-----------------------------------------------------------------------------------------
class EpiPred(object):
    """
    Main class
    """

    def __init__(self,PeptidePath,Alleles,Method,ReportFile):
        """
        Initial setting for the commandline parameters
        """
        self.PeptidePath = PeptidePath
        self.Alleles = Alleles
        self.Method = Method
        self.ReportFile = ReportFile
        self.FormatedPath = self.PeptidePath.rstrip('/')
    
    def ValidateFasta(self,FastaFile):
        """
        Determine peptide size and check for empty or improperly formated
        fasta files
        """
        SizeCount = 0
        f = open(FastaFile, "r")
        StartCharacters = dict()
        EmptyLine = re.compile('^$')
        for line in f:
            if not re.match(EmptyLine,line):
                StrippedLine = line.strip()
                StartCharacter = StrippedLine[0]
                StartCharacters[StartCharacter] = ''
                
                if not StrippedLine.startswith('>'):
                    LineSize = len(StrippedLine)
                    SizeCount += LineSize
        if '>' not in StartCharacters or SizeCount == 0:
            sys.exit("\nError!: Empty or improperly formated fasta file encountered\n")
        else:
            pass

    def Predict(self):
        """
        Run the prediction
        """
        Peptides = os.listdir(self.FormatedPath)
        for Allele in self.Alleles:
            for Peptide in Peptides:
                PeptidePath = self.FormatedPath + "/" + Peptide
                self.ValidateFasta(PeptidePath)
                PredResults = self.FormatedPath + "/" + Peptide + "_" + Allele + ".txt"
                if os.path.exists(PredResults):
                    pass
                else:
                    os.system("mhc_II_binding.py %s %s %s > %s" % \
                    (self.Method,Allele,PeptidePath,PredResults))


    
    def Parseconsensus(self,PredFile):
        """
        Parse consensus predicted results
        """
        FoundPreds = [];
        f = open(PredFile);
        for line in f:
            if not line.startswith("allele"):
                SplitLine = line.strip().split('\t')
                Al = SplitLine[0]
                SeqNum = SplitLine[1]
                Start = SplitLine[2]
                End = SplitLine[3]
                Pep = SplitLine[4]
                ConRank = float(SplitLine[5])
                PredHolder = dict()
                    #if ConRank <= 1.0:
                PredLine = SeqNum+':'+Pep+':'+str(ConRank)
                PredHolder[Al] = PredLine
                FoundPreds.append(PredHolder)

        return FoundPreds
    


    def GenReport(self):
        """
        Generate prediction report
        """
        Predictions = glob.glob(self.FormatedPath + '/*.txt')
        Strains = []
        Alleles = []
        MatRaw = []
        for p in Predictions:
            DataArray = self.Parseconsensus(p)
            Strain = p.split('/')[-1].split('_')[0]
            Allele = p.split('/')[-1].split('_')[1].split('.')[0]
            UniqPeptides = []
            for i in DataArray:
                for key,value in i.items():
                    UniqPeptides.append(value.split(':')[1])
            UniqPeptides = list(set(UniqPeptides))
            PeptideCount = len(UniqPeptides)

            Strains.append(Strain)
            Alleles.append(Allele)
            MatRaw.append((Strain,PeptideCount,UniqPeptides))

        UniqAlleles = list(set(Alleles))
        #Header details
        HeadStrain = "Strain"
        HeadAlleles = "\t".join([i+'_Allele' for i in UniqAlleles])
        HeadPeptides = "\t".join([i+'_Peptide' for i in UniqAlleles])

        FOut = open(self.ReportFile,'w')
        FOut.write("%s\t%s\t%s\n" % (HeadStrain,HeadAlleles,HeadPeptides))
        for s in list(set(Strains)):
            FoundCount = []
            FoundPep = []
            for m in MatRaw:
                st,pc,unip = m
                if st == s:
                    FoundCount.append(pc)
                    FoundPep.append(unip)

            FlatPep = []
            for i in FoundPep:
                FlatPep.append(",".join(i))

            FOut.write("%s\t%s\t%s\n" % (s,"\t".join(map(str,FoundCount)),"\t".join(FlatPep)))



#ListReader
#-----------------------------------------------------------------------------------------
def ListReader(infile):
    Alleles = []
    for line in open(infile,"rU"):
        Alleles.append(line.strip())
    return Alleles


#Commandline Options
#-----------------------------------------------------------------------------------------
def CommandlineOptions():
    
    Parser = OptionParser(usage="usage: %prog [options]")
    Parser.add_option("-p", "--peptides",
                      type="string",
                      dest="PeptidePath",
                      help="A path to the directory holding the peptides")
    Parser.add_option("-a", "--allele",
                    type="string",
                    dest="AlleleFile",
                    help="File with list of alleles")
    Parser.add_option("-m", "--method",
                    type="string",
                    dest="Method",
                    default="consensus",
                    help="Prediction method")

    Parser.add_option("-o", "--output",
                    action="store",
                    type="string",
                    dest="ReportFile",
                    default="report.txt",
                    help="The final output file")

    (Options, Args) = Parser.parse_args()
    OptionsArgsParser = [Options,Args,Parser]
    
    return OptionsArgsParser

#Dependency check
#-----------------------------------------------------------------------------------------
def Which(program):
    """
    Check if a program has been installed in a *NIX system and is in the path
    """
    status = 0
    try:
        null = open("/dev/null", "w")
        subprocess.Popen(program, stdout=null, stderr=null)
        null.close()
        status = 1
    
    except OSError:
        status = 0
    
    return status



#Core Prog
#-----------------------------------------------------------------------------------------
def main():
    Options, Args, Parser = CommandlineOptions()
    PeptidePath = Options.PeptidePath
    AlleleFile = Options.AlleleFile
    Method = Options.Method
    ReportFile = Options.ReportFile

    if AlleleFile == None or PeptidePath == None:
        sys.stderr.write("\nError: A mandatory option is missing!\n")
        Parser.print_help()
        sys.exit(-1)

    PredicStatus = Which('predict_binding.py')
    if PredicStatus == 0:
        sys.exit("\nError!: MHC predictor is not installed\n")
    else:
        Alleles = ListReader(AlleleFile)
        PredRunner = EpiPred(PeptidePath,Alleles,Method,ReportFile)
        PredRunner.Predict()
        PredRunner.GenReport()

#Run it
#-----------------------------------------------------------------------------------------
if __name__ == '__main__':
    main()