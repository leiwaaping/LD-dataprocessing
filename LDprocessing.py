#!/usr/bin/python3
# -*- coding: UTF-8 -*-
# @Date    : 2019-3-7
# @contact  : beccalihp@yahoo.com
# @Version : 2.0
# @input format: vcfs in one dir, target_snp list, distance
# @command :python *.py -t targetsnp.list -i vcfsDir -o output.file -d distance between loci

import sys, getopt, re
from typing import Any, Union
from collections import Counter
import os
import datetime
from sys import stdin
import linecache
import re

def main(argv):
   starttime = datetime.datetime.now()
   try:
      opts, args = getopt.getopt(sys.argv[1:],"ht:i:o:d:")
      inputfile = ''
      outputfile = ''
      targetfile = ''
      LDmaxLen = ''
   except getopt.GetoptError:
      print('An error has occured. \nYou can use this by: python compro.py -t <targetfile> -i <vcfdir> -o <outputfile> -d <maxdistance>')
      sys.exit(2)
   for opt, arg in opts:
      if opt == '-h':
         print('python *.py -t <targetfile> -i <inputfile> -o <outputfile> -o <distance>')
         sys.exit()
      elif opt == "-t":
         targetfile = arg
      elif opt == "-i":
         inputfile = arg
      elif opt == "-o":
         outputfile = arg
      elif opt == "-d":
         LDmaxLen = arg
   print("target snp list in :",targetfile)
   print("vcfs saved in dir :",inputfile)
   print("output file :",outputfile)
   print("max distance between loci is :",LDmaxLen)
   compro(targetfile,inputfile)
   cout_candidate_loci(LDmaxLen,outputfile)
   #show runtime
   endtime = datetime.datetime.now()
   time = "Run time is %s s\n"%(endtime - starttime).seconds
   print (time)

def file_name(file_dir):   
    L=[]   
    for dirpath, dirnames, filenames in os.walk(file_dir):  
        for file in filenames :  
            if os.path.splitext(file)[1] == '.vcf':  
                L.append(os.path.join(dirpath, file))  
    return L 

def compro(targetfile,inputfile):
    print("################### processing compro #####################")
    snps = openfile(targetfile)
    vcffiles=file_name(inputfile)
    vcfsall=[]
    out = 'rs\tchr\tposition\tref\t'
    for i in range(len(vcffiles)):
        vcfs = openfile(vcffiles[i])
        vcfsall.append(vcfs)
        out+='\t'+vcffiles[i] #readin first line and add filename
    out += '\tmain_allele'
    for snp in snps:
        snp = snp.rstrip('\n')
        snpstr=""
        alt=[]
        ref=""
        for vcfs in vcfsall:# readin each file
            for vcf in vcfs:# readin each line
                matchObj = re.match(r"(.*)\t(.*)\t" + snp + "\t(.*)\t([A-Z]*)\t", vcf)
                if matchObj:
                    snpstr = snp + '\t' + matchObj.group(1) + '\t' + matchObj.group(2) + '\t' + matchObj.group(3) + '\t'
                    ref = matchObj.group(3)+'\t'
        if snpstr != "":
            out += '\n' + snpstr
            for vcfs in vcfsall:
                vcfnum1 = ""
                for vcf in vcfs:
                    match1 = re.match(r"(.*)\t(.*)\t" + snp + "\t(.*)\t([A-Z]*)\t", vcf)
                    if match1:
                        vcfnum1 = match1.group(4) + '\t'
                        alt.append(vcfnum1)
                        out += vcfnum1
                if vcfnum1=="":#if nnot match
                    out += ref
        count1 = count(alt)#every rs# with one max
        out += count1
    fw = open('tempunit.txt', 'w', encoding='utf-8')
    fw.write(out)
    fw.close

def openfile(filename):
    f=open(filename, "r+",encoding='utf-8')
    strs=f.readlines()
    f.close
    return strs

def count(alt):
    #define an dict
    d ={}
    #record the most common
    max_key = ""
    for i in alt:
        if i not in d:
            count = alt.count(i)
            d[i] = count
            if count > d.get(max_key, 0):
                max_key = i
    return max_key
    #print(max_count(alt))

def cout_candidate_loci(LDmaxLen,outputfile):
    print("################### counting candidate loci #####################")
    fout = open(outputfile,'w')
    fout.write("SNP1 SNP2 #AB #A #a #B #b #Sample\n")
    #calculating candidate LD position
    f = open('tempunit.txt','r')
    next(f)
    #k is recording line number of main loop
    k=1
    for line in f:
      #print(line)
    
      list1 = line.strip().split("\t")  
      #len of array,number of sample number = len -5
      sampleindex=int(len(list1))
      #print(sampleindex)
      #str = "%s\n"%list1
      #print(str)  
      Chr1 = list1[1]
      pos1 = int(list1[2])
      rs1 = list1[0]
      #print(Chr1,pos1,rs1)
      #linecache.checkcache
      #print(list1)
      #markdown which line we are in
      k = k+1
      #gothrough next pos the calculate distance between two loci
      j = k+1
      while j:
        line2 = linecache.getline('tempunit.txt', j)   
        if(line2 !=''):
          list2 = line2.strip().split("\t")
          Chr2 = list2[1]
          pos2 = int(list2[2])
          rs2 = list2[0]
          d =pos2-pos1
        else:
          break   
        #different chromosome
        if((Chr1!=Chr2)or(d>int(LDmaxLen))):
          #line = f.readline()
          #str = "%s  %s  %s break!\n"%(rs1,rs2,d)
          #print(str)
          break
        #else d < LDminLen,candidate LD loci
        else:
          AB = 0
          mainbase1 = list1[3]
          minbase1 = list1[len(list1)-1]
          #calculate p(A)and p(a)
          A = re.compile(r"%s"%list1[3])
          nA = len(re.findall(A,line))-1
          a = re.compile(r"%s"%list1[len(list1)-1])
          na = len(re.findall(a,line))
          #print(nA,na)
        
          mainbase2 = list2[3]
          minbase2 = list2[len(list2)-1]
          #calculate p(B)and p(b)
        
          ####python
          B = re.compile(r"%s"%list2[3])
          nB = len(re.findall(B,line2))-1
          b = re.compile(r"%s"%list2[len(list1)-1])
          nb = len(re.findall(b,line2))
          #print(nB,nb)
            
        
          for x in range(5,sampleindex-1):
            genotype1=list1[x]
            genotype2=list2[x]
            #print(genotype)
          
            #python
            L13 = re.compile(r"%s"%list1[3])
            L23 = re.compile(r"%s"%list2[3])
            if((re.match(L13,genotype1))and(re.match(L23,genotype2))):
              AB = AB + 1
         
          
          #str = "%s  %s  %s\n"%(rs1,rs2,d)
          samplenum = sampleindex - 5
          str = "%s %s %s %s %s %s %s %s\n"%(rs1,rs2,AB,nA,na,nB,nb,samplenum)
          fout.write(str)
          #print(str)
          j=j+1
    
    f.close

if __name__ == "__main__":
   main(sys.argv[1:])
