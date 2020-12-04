# kturnSeeker
kturnSeeker: A computational software for identifying forward Kink-turn RNAs (fktRNAs) and backward Kink-turn RNAs (bktRNAs) from fasta file.

Overview:
---------
kturnSeeker is an software to identify novel forward Kink-turn RNAs (fktRNAs) and backward Kink-turn RNAs (bktRNAs) from fasta file.

Usage:
---------
Usage:  kturnSeeker [options] --fa <fasta file>
[options]
-o/--output <file>     : output file
-v/--verbose           : verbose information
-V/--version           : kturnSeeker version
-h/--help              : help informations
-i/--min-dist          : minimum distance between two kturn elements[default=15]
-a/--max-dist          : maximum distance between two kturn elements[default=500]
-m/--min-score         : minimum score for C-stem pairs of kturn[default>=5]
-l/--min-len           : minimum length for C-stem pairs of kturn[default>=2]
-f/--nofw              : will not scan forward reference strand
-r/--norc              : will not scan reverse-complement reference strand
-s/--stringent         : output strigent for stem pairs of kturn


Installation:<BR>
---------
Download kturnSeeker-0.1.tar.gz from https://github.com/sysu-software/kturnSeeker/releases ; unpack it, and make:<BR>
tar -xzvf kturnSeeker-0.1.tar.gz<BR>
cd kturnSeeker-0.1<BR>
make<BR>

System requirements:
---------
Operating system: kturnSeeker is designed to run on POSIX-compatible platforms, including UNIX, Linux and Mac OS/X. We have tested  most extensively on Linux and MacOS/X because these are the machines we develop on.<BR>
Compiler: The source code is compiled with  the C++ compiler g++. We test the code using the g++ compilers.<BR>
By default, kturnSeeker does not require any additional libraries to be installed by you.<BR>

Prerequisites:<BR>
---------
The input of kturnSeeker is a FASTA file<BR>

run kturnSeeker:
---------
cd example;
./kturnSeeker -m 7 --fa  SNORD126.fa \>kturnSeeker_candidates.txt<BR>

Output:
---------
#seqName	start	end	kturnName	score	strand	kturnType	stemLen	seqLen	5pDist	3pDist	kturnSeq	orgSeq	stem1Seq	pairs	stem2Seq
SNORD126	0	77	forward_kturn-1	15	+	forward	7	77	0	0	AGTTTGCCATGATGAAATGCATGTTAAGTCCGTGTTTCAGCTGATCAGCCTGATTAAACACATGCTCTGAGCAGACT	AGTTTGCCATGATGAAATGCATGTTAAGTCCGTGTTTCAGCTGATCAGCCTGATTAAACACATGCTCTGAGCAGACT	stem1: 5'-AGTTTGCCATGATG\	pairs:    |||:|||   xx.| AAATGCATGTTAAGTCCGTGTTTCAGCTGATCAGCCTGATTAAACACATGCT	stem2: 3'-TCAGACG---AGTC/

Note: # is comment line<BR>

Acknowledgements:
---------
Thanks a lot to everyone who contributed to the public used by kturnSeeker.<BR>

Contact :
---------
*****************************************************************************************<BR>
 \*	kturnSeeker - A computational software for identifying Kink-turn RNAs.<BR>
 \*<BR>
 \*	Author : Jian-Hua Yang <yangjh7@mail.sysu.edu.cn><BR>
 \* <BR>
 \*	RNA Information Center, School of Life Sciences, Sun Yat-Sen University<BR>
 \*	<BR>
 \*  Create date: 11/18/2019<BR>
 \*  <BR>
 \*  last modified time: 09/01/2020<BR>
 ****************************************************************************************<BR>

