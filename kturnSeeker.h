/********************************************************
 * kturnSeeker: identify kink-turn RNAs from fasta file
 * jianhua yang yangjh7@mail.sysu.edu.cn
 * $ 2019/9/18
 *******************************************************/
#ifndef kturnSeeker_HEAD_H
#define kturnSeeker_HEAD_H
#include "fasta.h"

#define gcScore 3
#define auScore 2
#define guScore -5
#define misScore -5

struct parameterInfo {
	int verbose;
	int nofw;
	int norc;
	int stringent;
	unsigned int maxDist;
	unsigned int minDist;
	double minStemScore;
	int minStemLen;
};

typedef struct parameterInfo parameter;

struct stemScoreInfo {
	int length;
	double score;
};

typedef struct stemScoreInfo stemScore;

void searchInFastaFile(FILE *outfp, parameter *paraInfo, FILE *fafp);

int kturnScan(FILE *outfp, parameter *paraInfo, FASTA *fastaSeq, char strand, int index);
int outputRevKturn(FILE *outfp, parameter *paraInfo, int start, int end, FASTA *fastaSeq, int externPairs, char strand, int index);

int outputForKturn(FILE *outfp, parameter *paraInfo, int start, int end, FASTA *fastaSeq, int externPairs, char strand, int index);

int overlap(int s1, int e1, int s2, int e2);

int drawRevPairs(char *seq1, char *seq2, char *pairs, int stemLen);

int drawForPairs(char *seq1, char *seq2, char *pairs, int fixLen);

int structureFilter(char *structure, int stemLen);

int getForCstemLen(int start, int end, char *seq, int seqLen, stemScore *ss);

int getRevCstemLen(int start, int end, char *seq, int seqLen, stemScore *ss);

int getGenomePos(FILE *outfp, char *name, int chromStart, int chromEnd, int index, int score, const char *type);

int encodeInt(char ch);

int RNApair(char bp1, char bp2);

int GApair(char bp1, char bp2);

#endif /* End kturnSeeker_HEAD_H */
