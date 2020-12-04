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

void searchInFastaFile(FILE *outfp, parameter *paraInfo, FILE *fafp)
{
	FASTA fastaSeq;
	int i = 0;
	char strand = '+';
	int kturnLen = 15;
	int index = 1;
	fprintf(outfp, "#seqName\tstart\tend\tkturnName\tscore\tstrand\tkturnType\tstemLen\tseqLen\t5pDist\t3pDist\tkturnSeq\torgSeq\tstem1Seq\tpairs\tstem2Seq\n");
	while (readSeq(fafp, &fastaSeq) != 0)
	{
		if (paraInfo->verbose) fprintf(stderr, "seqName: %s\n", fastaSeq.seqName);
		fflush(outfp);
		if (strlen(fastaSeq.seq) < paraInfo->minDist + 15)
		{
			freeFasta(&fastaSeq);
			continue;
		}
		if (strand == '+') {
			if (!paraInfo->nofw)
				index = kturnScan(outfp, paraInfo, &fastaSeq, strand, index);
			strand = '-';
		}
		if (strand == '-') {
			if (!paraInfo->norc)
			{
				reverseComp(fastaSeq.seq);
				index = kturnScan(outfp, paraInfo, &fastaSeq, strand, index);
			}
			strand = '+';
		}
		freeFasta(&fastaSeq);
		i++;
	}
	fprintf(stderr, "identify %d kturn structure\n", index - 1);
}

int kturnScan(FILE *outfp, parameter *paraInfo, FASTA *fastaSeq, char strand, int index)
{
	/* kturn structure

	forward structure
	5'-GCTGTTCGTTGATG\
	   |||||||   ...| AGGCTCAGAGTGAGCGCTGGGTACAGCGCCCGAATCGGACAGTGTAGAACCATTCTCTACTGCCTTCCTT
	3'-CGACAAG---AGTC/

	backward structure
	GCCTCTGA---GCCTCATCTATAGAAGGAGGCTATGATGAA
	    |...   |||||
	  AAGTAGTATCGGAG

	*/
	int i = 0;
	int j = 0;
	int bulgeLen = 3;
	int k        = 0;
	int minKturnLen = 18;
	int minStemLen  = 3;
	int kturnEleLen = 4; // length of second element
	int seqLen = fastaSeq->seqLen;
	char *seq  = fastaSeq->seq;
	int start = 0;
	int end   = 0;
	int secEleStart = 0; // start of secondary element
	int secEleEnd   = 0; // end of secondary element
	int bestStart   = -1000;
	int bestEnd     = -1000;
	int bestExStemLen = -1;
	int findKturnTag  = 0;
	int bestKturnType = 0;
	int kturnType     = 0;
	int exStemLen     = 0;
	double score      = 0;
	double bestStemScore = -5;
	stemScore *ss = (stemScore *)safeMalloc(sizeof(stemScore));
	ss->length = 0;
	ss->score  = 0;
	for (i = 0; i < seqLen - minKturnLen; i++)
	{
		if (toupper(seq[i + 2]) == 'G' // G.A
		        && toupper(seq[i + 3]) == 'A') // G.A
		{
			secEleStart = i + paraInfo->minDist;
			secEleEnd = secEleStart + paraInfo->maxDist;
			for (j = secEleStart; j < secEleEnd && j < seqLen - kturnEleLen; j++)
			{
				exStemLen    = 0;
				kturnType    = 0;
				ss->length   = 0;
				ss->score    = 0;
				score        = 0;

				if (toupper(seq[j]) == 'G' // G.A
				        && toupper(seq[j + 1]) == 'A' //G.A
				   )
				{
					// search for backward NC-stem
					if (
					    RNApair(seq[i], seq[j + 3]) >= 3 // NC-stem G-C (only) and A-U pairs, 2020/06/20
					    && RNApair(seq[i + 1], seq[j + 2]) <= 1 // G-U or mismatches
					    //&& toupper(seq[i + 0]) == 'C' // CTGA motif
					    //&& toupper(seq[i + 1]) == 'T' // CTGA motif
					)
					{
						start = i + 4;
						end   = j - 4;
						if (RNApair(seq[start], seq[end]) >= 2
						        && RNApair(seq[start + 1], seq[end - 1]) >= 2
						        && RNApair(seq[start + 2], seq[end - 2]) >= 1 // for U4ATAC others is G-C or A-U pairs
						   )
						{
							exStemLen = getRevCstemLen(start, end, seq, seqLen, ss);
							score = ss->score;
							if (exStemLen >= paraInfo->minStemLen && ss->score >= paraInfo->minStemScore)
							{
								kturnType = 2;
							}
						}
					}
					// search for forward NC-stem
					if (i >= 5 // need minim five bases for stem and bulge
					        && RNApair(seq[i + 5], seq[j - 2]) >= 3 // NC-stem G-C (only) and A-U pairs, 2020/06/20
					        && RNApair(seq[i + 4], seq[j - 1]) <= 1 // NC-stem, G-U or mismatches
					        //&& toupper(seq[j - 2]) == 'C' // CTGA motif
					        //&& toupper(seq[j - 1]) == 'T' // CTGA motif
					   )
					{
						start = i - 2;
						end   = j + 2;
						if (RNApair(seq[start], seq[end]) >= 2
						        && RNApair(seq[start - 1], seq[end + 1]) >= 2
						        && RNApair(seq[start - 2], seq[end + 2]) >= 1 // for U4ATAC others is G-C or A-U pairs
						   )
						{
							exStemLen = getForCstemLen(start, end, seq, seqLen, ss);
							score = ss->score;
							if (exStemLen >= paraInfo->minStemLen && ss->score >= paraInfo->minStemScore)
							{
								kturnType = 1;
							}
						}
					}
				}
				if (kturnType == 0) continue;
				if (overlap(start, end, bestStart, bestEnd) > 0)
				{
					if (score >= bestStemScore)
					{
						bestExStemLen = exStemLen;
						bestStart     = start;
						bestEnd       = end;
						bestKturnType = kturnType;
						bestStemScore = score;
					}
				}
				else
				{
					if (bestExStemLen >= paraInfo->minStemLen && bestStemScore >= paraInfo->minStemScore)
					{
						if (paraInfo->verbose) fprintf(stderr, "output kturn %s StemLen:%d StemScore:%.5f\n", fastaSeq->seqName, bestExStemLen, bestStemScore);
						if ( bestKturnType == 1)
							index = outputForKturn(outfp, paraInfo, bestStart, bestEnd, fastaSeq, bestExStemLen, strand, index);
						if ( bestKturnType == 2)
							index = outputRevKturn(outfp, paraInfo, bestStart, bestEnd, fastaSeq, bestExStemLen, strand, index);
					}
					bestExStemLen = exStemLen;
					bestStart = start;
					bestEnd = end;
					bestKturnType = kturnType;
					bestStemScore = score;
				} //
			} // for j
		} // if G.A
	}// for i
	// last time output
	if (bestExStemLen >= paraInfo->minStemLen && bestStemScore >= paraInfo->minStemScore)
	{
		if (paraInfo->verbose) fprintf(stderr, "output kturn %s StemLen:%d StemScore:%.5f\n", fastaSeq->seqName, bestExStemLen, bestStemScore);
		if ( bestKturnType == 1)
			index = outputForKturn(outfp, paraInfo, bestStart, bestEnd, fastaSeq, bestExStemLen, strand, index);
		if ( bestKturnType == 2)
			index = outputRevKturn(outfp, paraInfo, bestStart, bestEnd, fastaSeq, bestExStemLen, strand, index);
	}
	safeFree(ss);
	return index;
}

int getRevCstemLen(int start, int end, char *seq, int seqLen, stemScore *ss)
{
	int i = 0;
	int externPairs = 0;
	int maxExternLen = 15;
	int bestStemScore = 0;
	int bestExStemLen = 0;
	int stemScore = 0;
	int externLen = end;
	for (i = 0; i < maxExternLen && i < externLen; i++)
	{
		if (start + i >= end - i) break;
		int pair = RNApair(seq[start + i], seq[end - i]);
		if (pair)
		{
			if (pair == 1) stemScore += guScore;
			if (pair == 2) stemScore += auScore;
			if (pair == 3) stemScore += gcScore;
			if (stemScore > bestStemScore)
			{
				bestStemScore = stemScore;
				externPairs = i + 1;
			}
		}
		else
		{
			stemScore += misScore;
		}
	}
	ss->score = bestStemScore;
	ss->length = externPairs;
	return externPairs;
}

int getForCstemLen(int start, int end, char *seq, int seqLen, stemScore *ss)
{
	int i = 0;
	int externLen = start;
	int externPairs = 0;
	int maxExternLen = 15;
	int bestStemScore = 0;
	int bestExStemLen = 0;
	int stemScore = 0;
	if (externLen > seqLen - end) externLen = seqLen - end;
	for (i = 0; i <= externLen && i < maxExternLen; i++)
	{
		int pair = RNApair(seq[start - i], seq[end + i]);
		if (pair)
		{
			if (pair == 1) stemScore += guScore;
			if (pair == 2) stemScore += auScore;
			if (pair == 3) stemScore += gcScore;
			if (stemScore > bestStemScore)
			{
				bestStemScore = stemScore;
				externPairs = i + 1;
			}
		}
		else
		{
			stemScore += misScore;
		}
	}
	ss->score = bestStemScore;
	ss->length = externPairs;
	return externPairs;
}

int overlap(int s1, int e1, int s2, int e2)
// overlapping bases
{
	int maxMin = s1;
	int minEnd = e1;
	if (s2 > maxMin) maxMin = s2;
	if (e2 < minEnd) minEnd = e2;
	return minEnd - maxMin;
}

int outputRevKturn(FILE *outfp, parameter *paraInfo, int start, int end, FASTA *fastaSeq, int externPairs, char strand, int index)
{
	/*
	GCCTCTGA---GCCTCATCTATAGAAGGAGGCTATGATGAA
	    |...   |||||
	  AAGTAGTATCGGAG
	start = i + 4;
	end   = j - 4;
	*/
	int fixStemLen = 4;
	int stemBulge = 3;
	char *seq = fastaSeq->seq;
	start = start - fixStemLen;
	end = end + fixStemLen + stemBulge + 1;
	int seqLen = end - start;
	int stemLen = externPairs;
	int kturnLen = stemLen + fixStemLen + stemBulge;
	int stem1Len = stemLen + fixStemLen;
	int stem2Len = stemLen + fixStemLen + stemBulge;
	int loopLen = seqLen - stem1Len - stem2Len;
	int score = 0;
	int chromStart = start;
	int chromEnd = end;
	double mfe = 0;
	char stemBulgeSeq[] = "---";
	char *stem1Seq = (char *)safeMalloc(sizeof(char) * (kturnLen + 1));
	char *stem2Seq = (char *)safeMalloc(sizeof(char) * (kturnLen + 1));
	char *stem3Seq = (char *)safeMalloc(sizeof(char) * (stemLen + 1));
	char *pairStr = (char *)safeMalloc(sizeof(char) * (kturnLen + 1));
	char *loopSeq = (char *)safeMalloc(sizeof(char) * (loopLen + 1));
	char *structure = (char *)safeMalloc(sizeof(char) * (seqLen + 1));
	char repChar;
	strncpy(stem1Seq, seq + start, fixStemLen);
	stem1Seq[fixStemLen] = '\0';
	strcat(stem1Seq, stemBulgeSeq);
	strncpy(stem3Seq, seq + start + fixStemLen, stemLen);
	stem3Seq[stemLen] = '\0';
	strcat(stem1Seq, stem3Seq);

	strncpy(stem2Seq, seq + end - kturnLen, kturnLen);
	stem2Seq[kturnLen] = '\0';
	reverseBytes(stem2Seq, kturnLen);

	strncpy(loopSeq, seq + start + stem1Len, loopLen);
	loopSeq[loopLen] = '\0';
	score = drawForPairs(stem1Seq, stem2Seq, pairStr, fixStemLen + stemBulge);
	repChar = seq[end];
	seq[end] = '\0';
	int tag = 0;
	//mfe = fold(seq + start, structure);
	if (score >= paraInfo->minStemScore)
		//&& structureFilter(structure, stemLen))
	{
		if (strand == '-')
		{
			chromStart = fastaSeq->seqLen - chromEnd;
			chromEnd = chromStart + seqLen;
		}
		//fprintf(stderr, "%s\n", fastaSeq->seqName);
		tag = getGenomePos(outfp, fastaSeq->seqName, chromStart, chromEnd, index, score, "backward");
		if (tag == 0)
			fprintf(outfp, "%s\t%d\t%d\tbackward_kturn-%d\t%d\t%c\tbackward\t", fastaSeq->seqName, chromStart, chromEnd, index, score, strand);
		fprintf(outfp, "%d\t%ld\t%d\t%ld\t", stemLen, fastaSeq->seqLen, chromStart, fastaSeq->seqLen - chromEnd);
		fprintf(outfp, "%s\t", seq + start);
		seq[end] = repChar;
		fprintf(outfp, "%s\t", seq);
		//fprintf(outfp, "%s\t", structure);
		// output stem pairs
		fprintf(outfp, "stem1: 5'-%s\\\t", stem1Seq);
		fprintf(outfp, "pairs:    %s %s\t", pairStr, loopSeq);
		fprintf(outfp, "stem2: 3'-%s/\n", stem2Seq);
		index++;
	}
	seq[end] = repChar;
	fflush(outfp);
	safeFree(stem1Seq);
	safeFree(stem2Seq);
	safeFree(loopSeq);
	safeFree(pairStr);
	safeFree(stem3Seq);
	safeFree(structure);
	return index;
}

int outputForKturn(FILE *outfp, parameter *paraInfo, int start, int end, FASTA *fastaSeq, int externPairs, char strand, int index)
{
	/* kturn structure
	5'-GCTGTTCGTTGATG\
	   |||||||   ...| AGGCTCAGAGTGAGCGCTGGGTACAGCGCCCGAATCGGACAGTGTAGAACCATTCTCTACTGCCTTCCTT
	3'-CGACAAG---AGTC/
	start = i - 2;
	end   = j + 2;
	*/
	char *seq = fastaSeq->seq;
	start = start - externPairs + 1;
	end = end + externPairs;
	int seqLen = end - start;
	int fixStemLen = 4;
	int stemBulge = 3;
	int stemLen = externPairs;
	int kturnLen = stemLen + fixStemLen + stemBulge;
	int stem1Len = stemLen + fixStemLen + stemBulge;
	int stem2Len = stemLen + fixStemLen;
	int loopLen = seqLen - stem1Len - stem2Len;
	int score = 0;
	int chromStart = start;
	int chromEnd = end;
	double mfe = 0;
	char stemBulgeSeq[] = "---";
	char *stem1Seq = (char *)safeMalloc(sizeof(char) * (kturnLen + 1));
	char *stem2Seq = (char *)safeMalloc(sizeof(char) * (kturnLen + 1));
	char *stem3Seq = (char *)safeMalloc(sizeof(char) * (stemLen + 1));
	char *pairStr = (char *)safeMalloc(sizeof(char) * (kturnLen + 1));
	char *loopSeq = (char *)safeMalloc(sizeof(char) * (loopLen + 1));
	char *structure = (char *)safeMalloc(sizeof(char) * (seqLen + 1));
	char repChar;
	strncpy(stem1Seq, seq + start, kturnLen);
	stem1Seq[kturnLen] = '\0';
	strncpy(stem2Seq, seq + end - stem2Len, fixStemLen);
	stem2Seq[fixStemLen] = '\0';
	strcat(stem2Seq, stemBulgeSeq);
	strncpy(stem3Seq, seq + end - stemLen, stemLen);
	stem3Seq[stemLen] = '\0';
	strcat(stem2Seq, stem3Seq);
	reverseBytes(stem2Seq, kturnLen);
	strncpy(loopSeq, seq + start + stem1Len, loopLen);
	loopSeq[loopLen] = '\0';
	score = drawRevPairs(stem1Seq, stem2Seq, pairStr, stemLen);
	repChar = seq[end];
	seq[end] = '\0';
	int tag = 0;
	//mfe = fold(seq + start, structure);
	if (score >= paraInfo->minStemScore)
		//&& structureFilter(structure, stemLen))
	{
		if (strand == '-')
		{
			chromStart = fastaSeq->seqLen - chromEnd;
			chromEnd = chromStart + seqLen;
		}
		//fprintf(stderr, "%s\n", fastaSeq->seqName);
		tag = getGenomePos(outfp, fastaSeq->seqName, chromStart, chromEnd, index, score, "forward");
		if (tag == 0)
			fprintf(outfp, "%s\t%d\t%d\tforward_kturn-%d\t%d\t%c\tforward\t", fastaSeq->seqName, chromStart, chromEnd, index, score, strand);
		fprintf(outfp, "%d\t%ld\t%d\t%ld\t", stemLen, fastaSeq->seqLen, chromStart, fastaSeq->seqLen - chromEnd);
		fprintf(outfp, "%s\t", seq + start);
		seq[end] = repChar;
		fprintf(outfp, "%s\t", seq);
		//fprintf(outfp, "%s\t", structure);
		// output stem pairs
		fprintf(outfp, "stem1: 5'-%s\\\t", stem1Seq);
		fprintf(outfp, "pairs:    %s %s\t", pairStr, loopSeq);
		fprintf(outfp, "stem2: 3'-%s/\n", stem2Seq);
		index++;
	}
	seq[end] = repChar;
	fflush(outfp);
	safeFree(stem1Seq);
	safeFree(stem2Seq);
	safeFree(loopSeq);
	safeFree(pairStr);
	safeFree(stem3Seq);
	safeFree(structure);
	return index;
}

// >chr1:15560083-15560125(+)
int getGenomePos(FILE *outfp, char *name, int chromStart, int chromEnd, int index, int score, const char *type)
{
	int tag = 0;
	char *chr = strchr(name, '-');
	if (chr != NULL)
	{
		int pos1 = chr - name;
		name[pos1] = ':';
		chr = strchr(name, '(');
		if (chr == NULL) return tag;
		int pos2 = chr - name;
		name[pos2] = ':';

		int fieldNum = 0;
		char **fields = NULL;
		char delims[] = ":";
		char strand = '+';
		int start = 0;
		int end = 0;
		char *targetChrom = NULL;
		int targetStart = 0;
		int targetEnd = 0;
		fields = splitString(name, delims, &fieldNum);
		if (fieldNum >= 4) {
			tag = 1;
			targetChrom = fields[0];
			start = atoi(fields[1]);
			end = atoi(fields[2]);
			strand = fields[3][0];
			if (strand == '+') {
				targetStart = start + chromStart;
				targetEnd = start + chromEnd;
			}
			else {
				targetStart = end - chromEnd;
				targetEnd = end - chromStart;
			}
			if (targetStart < 0) targetStart = 0;
			name[pos1] = '-';
			name[pos2] = '(';
			fprintf(outfp, "%s\t%d\t%d\t%s_kturn-%d:%s\t%d\t%c\t%s\t", targetChrom, targetStart, targetEnd, type, index, name, score, strand, type);
		}
		freeWords(fields, fieldNum);
		name[pos1] = '-';
		name[pos2] = '(';
	}
	return tag;
}

int structureFilter(char *structure, int stemLen)
{
	int i = 0;
	int stemTag = 1;
	//short int *pairTable = make_pair_table(structure);
	//int structLen = strlen(structure);
	//for (i = 2; i <= stemLen && i < structLen; i++)
	//{
	//	if (pairTable[i] != (structLen - i + 1))
	//	{
	//		//fprintf(stderr, "%d %d\n", pairTable[i], structLen);
	//		stemTag = 0;
	//		break;
	//	}
	//}
	//safeFree(pairTable);
	return stemTag;
}

int drawForPairs(char *seq1, char *seq2, char *pairs, int fixLen)
{
	int i = 0;
	int minLen = strlen(seq1);
	int seq2Len = strlen(seq2);
	int score = 0;
	if (minLen > seq2Len) minLen = seq2Len;
	for (i = 0; i < minLen; i++)
	{
		int pair = RNApair(seq1[i], seq2[i]);
		if (pair > 0)
		{
			if (pair == 1)
			{
				pairs[i] = ':';
				if (i >= fixLen) score += 0;
			}
			if (pair == 2)
			{
				pairs[i] = '|';
				if (i >= fixLen) score += auScore;
			}
			if (pair == 3)
			{
				pairs[i] = '|';
				if (i >= fixLen) score += gcScore;
			}
		}
		else
		{
			pairs[i] = '.';
			if (seq1[i] == '-' || seq2[i] == '-') pairs[i] = ' ';
			if (GApair(seq1[i], seq2[i])) pairs[i] = 'x';
			if (i >= fixLen) score += misScore;
		}
	}
	return score;
}

int drawRevPairs(char *seq1, char *seq2, char *pairs, int stemLen)
{
	int i = 0;
	int minLen = strlen(seq1);
	int seq2Len = strlen(seq2);
	int score = 0;
	if (minLen > seq2Len) minLen = seq2Len;
	for (i = 0; i < minLen; i++)
	{
		int pair = RNApair(seq1[i], seq2[i]);
		if (pair > 0)
		{
			if (pair == 1)
			{
				pairs[i] = ':';
				if (i < stemLen) score += 0;
			}
			if (pair == 2)
			{
				pairs[i] = '|';
				if (i < stemLen) score += auScore;
			}
			if (pair == 3)
			{
				pairs[i] = '|';
				if (i < stemLen) score += gcScore;
			}
		}
		else
		{
			pairs[i] = '.';
			if (seq1[i] == '-' || seq2[i] == '-') pairs[i] = ' ';
			if (GApair(seq1[i], seq2[i])) pairs[i] = 'x';
			if (i < stemLen) score += misScore;
		}
	}
	return score;
}

int encodeInt(char ch)
{
	/* translate character to number */

	ch = toupper(ch);

	if (ch == 'A') return 0;
	if (ch == 'C') return 1;
	if (ch == 'G') return 2;
	if (ch == 'U' || ch == 'T') return 3;

	return 4;
}

int RNApair(char bp1, char bp2)
{
	int pairMatrix[5][5] = {
		/* A  C  G  T  N*/
		{0, 0, 0, 2, 0}, // A
		{0, 0, 3, 0, 0}, // C
		{0, 3, 0, 1, 0}, // G
		{2, 0, 1, 0, 0}, // T
		{0, 0, 0, 0, 0}  // N
	};
	return pairMatrix[(int)(encodeInt(bp1))][(int)(encodeInt(bp2))];
}

int GApair(char bp1, char bp2)
{
	int pairMatrix[5][5] = {
		/* A  C  G  T  N*/
		{0, 0, 1, 0, 0}, // A
		{0, 0, 0, 0, 0}, // C
		{1, 0, 0, 0, 0}, // G
		{0, 0, 0, 0, 0}, // T
		{0, 0, 0, 0, 0}  // N
	};
	return pairMatrix[(int)(encodeInt(bp1))][(int)(encodeInt(bp2))];
}
