# coding=utf-8
import time
import os
import sys
import subprocess

# read="/home/lab/gll/data_aln/bowtie2-master/example-hg19/reads/one.fastq"
# out1="/home/lab/gll/Effaln/example-hg19/index-32-r128/one.sam"

read="/home/lab/gll/data_aln/bowtie2-master/example-hg19/reads/ERR194146-1M.fastq"                         #Illumina 定长
out1="/home/lab/gll/Effaln-20220914/example-hg19/index-32-r128/1M.sam"
# read="/media/dell198/4fe135dd-6ee8-4bea-8ccd-f450dc5a7100/home/lab/xiaofei/data/SRR003177.recal.1M.fastq"    #LS454 变长
# out1="/home/lab/gll/mason_simulator/effaln-1M-ls454.sam"


exe1="/home/dell198/gll/Effaln-main-v1/effaln"
# exe1="/home/lab/gll/Effaln-20231011-lens/effaln"
# idx1="/home/lab/gll/Effaln-20231011-lens/example-hg19/index-32-r128/hg19"

def nnc(str1):
    p=subprocess.Popen(str1,shell=True)
    print(str1)
    p.communicate()

def test1():
    cmd=exe1+" -x "+idx1+" -U "+read+" -S "+out1
    # cmd=exe1+" -x "+idx1+" -f -U "+read+" -S "+out1
    start=time.time()
    nnc(cmd)
    end=time.time()
    take1=end-start
    print('running time: %f s' %(take1))
    return take1

if __name__=="__main__":
    # print('\n')
    print("----------testing 1  "+time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time())))
    take1=test1()