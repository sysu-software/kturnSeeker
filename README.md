# kturnSeeker
kturnSeeker: A computational software for identifying forward Kink-turn RNAs (fktRNAs) and backward Kink-turn RNAs (bktRNAs) from fasta file.

Overview:
---------
kturnSeeker is a software for identifying forward Kink-turn RNAs (fktRNAs) and backward Kink-turn RNAs (bktRNAs) from fasta file identified from RIP-PEN-seq/PEN-seq/sub-PEN-seq or others small RNA-seq data.

Usage:
---------
Usage:  kturnSeeker [options] --fa \<fasta file\><BR>
[options]<BR>
-o/--output <file>     : output file<BR>
-v/--verbose           : verbose information<BR>
-V/--version           : kturnSeeker version<BR>
-h/--help              : help informations<BR>
-i/--min-dist          : minimum distance between two kturn elements[default=15]<BR>
-a/--max-dist          : maximum distance between two kturn elements[default=500]<BR>
-m/--min-score         : minimum score for C-stem pairs of kturn[default>=5]<BR>
-l/--min-len           : minimum length for C-stem pairs of kturn[default>=2]<BR>
-f/--nofw              : will not scan forward reference strand<BR>
-r/--norc              : will not scan reverse-complement reference strand<BR>
-s/--stringent         : output strigent for stem pairs of kturn<BR>


Installation:<BR>
---------
Download kturnSeeker-1.0.tar.gz from https://github.com/sysu-software/kturnSeeker/releases ; unpack it, and make:<BR>
tar -xzvf kturnSeeker-1.0.tar.gz<BR>
cd kturnSeeker-1.0<BR>
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
../kturnSeeker -m 7 --fa  SNORD126.fa \>kturnSeeker_candidates.txt<BR>

Output:
---------
#seqName	start	end	kturnName	score	strand	kturnType	stemLen	seqLen	5pDist	3pDist	kturnSeq	orgSeq	stem1Seq	pairs	stem2Seq <BR>
SNORD126	0	77	forward_kturn-1	15	+	forward	7	77	0	0	AGTTTGCCATGATGAAATGCATGTTAAGTCCGTGTTTCAGCTGATCAGCCTGATTAAACACATGCTCTGAGCAGACT	AGTTTGCCATGATGAAATGCATGTTAAGTCCGTGTTTCAGCTGATCAGCCTGATTAAACACATGCTCTGAGCAGACT	stem1: 5'-AGTTTGCCATGATG\	pairs:    |||:|||   xx.| AAATGCATGTTAAGTCCGTGTTTCAGCTGATCAGCCTGATTAAACACATGCT	stem2: 3'-TCAGACG---AGTC/ <BR>

<STRONG>Note: # is comment line and annotation line</STRONG><BR>
<pre>
seqName\tSequence name
start\tStart position of K-turn structure in sequence.
end\tEnd position of K-turn structure in sequence.
kturnName	  Unique ID identification.
score       C-stem score of K-turn structure.
strand	     The direction of the strand.
kturnType   K-turn type (forward or backward).
stemLen	    The length of C-stem in K-trun structure.
seqLen      The length of sequence.
5pDist      The Distance between 5' motif and RNA 5' end.
3pDist	     The Distance between 3' motif and RNA 3' end.
orgSeq    	 The original RNA sequence.
kturnSeq	   The K-turn sequnce contained K-turn structure at ends of RNA.
stem1Seq	   The sequence of stem1(C-stem or NC-stem) in K-turn.
pairs	      The base-pairing in K-turn structure and loop sequences.
stem2Seq	   The sequence of stem2 (C-stem or NC-stem) in K-turn.
</pre>
<BR>
<STRONG>Note: You can paste these three columns (stem1Seq, pairs and stem2Seq) in three different rows. And then you can find that they are composed of a Kink-turn structure. For example:</STRONG><br>
forward kink-turn RNA<br>
<pre>
stem1: 5'-AGTTTGCCATGATG\
pairs:    |||:|||   xx.| AAATGCATGTTAAGTCCGTGTTTCAGCTGATCAGCCTGATTAAACACATGCT
stem2: 3'-TCAGACG---AGTC/
</pre>

backward kink-turn RNA<br>
<pre>
stem1: 5'-CTGA---GCACAG\	
pairs:    |.xx   |||||| AG	
stem2: 3'-GTAGTGACGTGTC/
</pre>
 
Acknowledgements:
---------
Thanks a lot to everyone who contributed to the public codes used by kturnSeeker.<BR>

Contact :
---------
*****************************************************************************************<BR>
 \*	kturnSeeker - A computational software for identifying Kink-turn RNAs.<BR>
 \*<BR>
 \*	Author : Jianhua Yang <yangjh7@mail.sysu.edu.cn><BR>
 \* <BR>
 \*	RNA Information Center, School of Life Sciences, Sun Yat-Sen University<BR>
 \*	<BR>
 \*  Create date: 11/18/2019<BR>
 \*  <BR>
 \*  last modified time: 04/10/2023<BR>
 ****************************************************************************************<BR>

