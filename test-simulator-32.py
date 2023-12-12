# 2023/12/8修改
# 测试不同长度的read
# 补充different error profile
# ablation analysis, 关于seed pruning

import time
import os
import sys
import subprocess

# len=sys.argv[1]                                                                  #read长度作为参数传入

# read="/home/lab/gll/simPan/mismatch_snp_reads_100M/sim_ms_1.fa.2wgsim"             # 1.MS
# out0="/home/dell198/gll/result-error-profile/effaln-ms.sam"
# samacc="/home/lab/gll/protest/mycode/MS-error-profile"

# read="/home/lab/gll/simPan/mismatch_reads_100M/sim_m_1.fa.2wgsim"             # 2.M
# out0="/home/dell198/gll/result-error-profile/effaln-m.sam"
# samacc="/home/lab/gll/protest/mycode/M-error-profile"

# read="/home/lab/gll/simPan/snp_reads_100M/sim_s_1.fa.2wgsim"             # 3.S
# out0="/home/dell198/gll/result-error-profile/effaln-s.sam"
# samacc="/home/lab/gll/protest/mycode/S-error-profile"

read="/home/lab/gll/simPan/reads_100M/sim_1.fa.2wgsim"             # 4.Exact
out0="/home/dell198/gll/result-error-profile/effaln-exact.sam"
samacc="/home/lab/gll/protest/mycode/Exact-error-profile"

# read="/home/lab/gll/simPan/mismatch_snp_reads_100M/sim_ms"                       #根据名称改
# out0="/home/dell198/gll/hisat2_simulator/effaln-1M-"+len+"-32.sam"
# samacc="/home/lab/gll/protest/mycode/samStatics_gold-1M-"+len                    #原本是针对mason2，这里要改为hisat2模拟器

# exe0="/home/lab/gll/Effaln-20231011-lens/effaln"
exe0="/home/dell198/gll/Effaln-main-v1/effaln"
idx0="/home/dell198/gll/calling/effaln.grch38"
# idx0="/home/lab/gll/Effaln-20231011-lens/example-hg19/index-32-r128/hg19"

def nnc(str1):
    p=subprocess.Popen(str1,shell=True)
    print(str1)
    p.communicate()

def test0():
    # cmd=exe0+" -x "+idx0+" -f -U "+read+" -S "+out0                    #调用
    cmd=exe0+" -x "+idx0+" -U "+read+" -S "+out0                    #调用
    start=time.time()
    nnc(cmd)
    end=time.time()
    take0=end-start
    print('running time: ########## %f s' %(take0))
    print(' ')
    nnc("/home/lab/gll/protest/mycode/samStaticsFlag4 "+out0)          #1.比对率
    print(' ')
    nnc(samacc+" "+out0)                                               #2.准确率
    return take0

if __name__=="__main__":
    # print('\n')
    print("----------testing 0  "+time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time())))
    take0=test0()
