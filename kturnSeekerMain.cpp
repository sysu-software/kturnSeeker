/********************************************************
 * kturnSeeker: identify kink-turn RNAs from fasta file
 * jianhua yang yangjh7@mail.sysu.edu.cn
 * $ 2019/9/18
 *******************************************************/
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<ctype.h>
#include<math.h>
#include <getopt.h>
#include <map>
#include <algorithm>
#include <ios>
#include <iostream>
#include <string>
#include <ostream>
#include <fstream>
#include <iomanip>
#include <locale>
#include <sstream>
#include <vector>
using namespace std;

#include "bioUtils.h"
#include "fasta.h"
#include "kturnSeeker.h"

void usage(void);
char version[] = "kturnSeeker version 0.1 [2015/05/12]";

int main(int argc, char **argv) {
  int i = 0;
  int c = 0;
  char *faFile  = NULL;
  FILE *fafp    = NULL;
  FILE *outfp   = NULL;
  char *outfile = NULL;
  int showVersion = 0;
  int showHelp = 0;
  parameter paraInfo;

  if (argc == 1) {
    usage();
  }

  const char *shortOptions = "vhfrsVo:i:a:m:F:l:";

  const struct option longOptions[] = {
    { "verbose" , no_argument , NULL, 'v' },
    { "help" , no_argument , NULL, 'h' },
    { "version" , no_argument , NULL, 'V' },
    { "fa" , required_argument , NULL, 'F' },
    { "output" , required_argument , NULL, 'o' },
    { "min-dist" , required_argument , NULL, 'i' },
    { "max-dist" , required_argument , NULL, 'a' },
    { "min-score" , required_argument , NULL, 'm' },
    { "min-len" , required_argument , NULL, 'l' },
    { "nofw" , no_argument , NULL, 'f' },
    { "norc" , no_argument , NULL, 'r' },
    { "stringent" , no_argument , NULL, 's' },
    {NULL, 0, NULL, 0} ,  /* Required at end of array. */
  };

  paraInfo.verbose = 0;
  paraInfo.minDist = 12;
  paraInfo.maxDist = 500;
  paraInfo.minStemScore = 5;
  paraInfo.minStemLen  = 2;
  paraInfo.nofw = 0;
  paraInfo.norc = 0;
  paraInfo.stringent = 0;

  while ((c = getopt_long(argc, argv, shortOptions, longOptions, NULL)) >= 0) {
    switch (c) {
    case 'v':
      paraInfo.verbose = 1;
      break;
    case 'h':
      showHelp = 1;
      break;
    case 'V':
      showVersion = 1;
      break;
    case 'f':
      paraInfo.nofw = 1;
      break;
    case 'r':
      paraInfo.norc = 1;
      break;
    case 's':
      paraInfo.stringent = 1;
      break;
    case 'o':
      outfile = optarg;
      break;
    case 'F':
      faFile  = optarg;
      break;
    case 'i':
      paraInfo.minDist = atoi(optarg);
      break;
    case 'a':
      paraInfo.maxDist = atoi(optarg);
      break;
    case 'm':
      paraInfo.minStemScore = atof(optarg);
      break;
    case 'l':
      paraInfo.minStemLen = atof(optarg);
      break;
    case '?':
      showHelp = 1;
      break;
    default:
      usage();
    }
  }

  if (faFile == NULL)
  {
    fprintf(stderr, "Please set the --fa option\n");
    exit(1);
  }
  fafp = (FILE *) fopen(faFile, "r");
  if (fafp == NULL) {
    fprintf(stderr, "ERROR: Can't open %s\n", faFile);
    fprintf(stderr, "Please set the --fa option\n");
    exit(1);
  }

  if (outfile == NULL) {
    outfp = stdout;
  } else {
    outfp = (FILE *) fopen(outfile, "w");
    if (outfp == NULL) {
      fprintf(stderr, "ERROR: Can't open the output file %s\n", outfile);
      usage();
    }
  }

  if (showVersion) {
    fprintf(stderr, "%s", version);
    exit(1);
  }

  if (showHelp) {
    usage();
    exit(1);
  }

  fprintf(stderr, "scan kturn motifs...\n");
  searchInFastaFile(outfp, &paraInfo, fafp);

  fflush(outfp);
  fclose(fafp);

  return 0;
}

void usage(void) {
  fprintf(stderr, "%s",
          "Usage:  kturnSeeker [options] --fa <fasta file>\n\
[options]\n\
-o/--output <file>     : output file\n\
-v/--verbose           : verbose information\n\
-V/--version           : kturnSeeker version\n\
-h/--help              : help informations\n\
-i/--min-dist          : minimum distance between two kturn elements[default=15]\n\
-a/--max-dist          : maximum distance between two kturn elements[default=500]\n\
-m/--min-score         : minimum score for C-stem pairs of kturn[default>=5]\n\
-l/--min-len           : minimum length for C-stem pairs of kturn[default>=2]\n\
-f/--nofw              : will not scan forward reference strand\n\
-r/--norc              : will not scan reverse-complement reference strand\n\
-s/--stringent         : output strigent for stem pairs of kturn\n\
");
  exit(1);
}
