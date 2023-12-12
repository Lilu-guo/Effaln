# Effaln
version 0.1 (<u>20220915</u>)   
version 0.2 (20231022)   
version 0.3 (20231205)

# What is it?
Effaln is an efficient Burrows-Wheeler-based mapper for longer Next-generation sequencing reads. 

# How to use it?
Effaln consists of two components, index building and read mapping. You should first build the FM-index with the reference genome, and then perform the mapping process.

## Step I. Install
  1. Download (or clone) the source code form https://github.com/guolilu3/Effaln
  2. Compile the source code. (Note that you need to compile semiWFA first)
     ```shell
     cd ./semiWFA   
     make
     
     cd ../
     make
     ```

## Step II. Build FM-index
  1. Run the shell command:
     ```shell
     ./effaln-index <refName> <idxName>
     ```
     where *refName* is the reference genome, *idxName* is the index file name.
  
## Step III. Mapping
  1. Run the shell command:
     ```shell
     ./effaln -x <idxName> -U <rdsName> -S <samName>
     ```
     where *idxName* is the index file name, *rdsName* is the sequencing reads file name, and *samName* is the mapping result file name.
  
# Feedback
Please report bugs to Email: guolilu@stu.xidian.edu.cn or jgxygll@163.com if any questions or suggestions.   
Your feedback and test results are welcome.

# License   
Effaln is available under the MIT license.   
