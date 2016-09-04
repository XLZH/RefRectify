/* This program is used to rectify the reference (.fa)
 * Input: the result of mpileup (.mp) and reference to be rectified (ref.fa)
 * Output: the corrected reference (out.fa)
 * Compile: gcc -std=c99 RefRectify.3.c -o RefRectify
 *
 * Version : 3
 * Time: 2016-8-21
 * Author: XLZH < xiaolongzhang2015@163.com >
 */

#define __STDC_WANT_LIB_EXT1__ 1
#include <stdio.h>
#include <linux/limits.h>
#include <getopt.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#define MAX_DEPTH 100000
#define MEM_UNIT 10000
#define TRANS(x, y) ( (x==0x2e || x==0x2c) ? toupper(y) : toupper(x) )

char Buff[MAX_DEPTH] = {'\0'};

typedef struct {      // eg: -12AAAAAAAAAAAA
    int digit;        // digit = 12
    int diglen;       // diglen = 2 (the length of string "12")
} DL;                 // Digit length

typedef struct {      // eg: T-5ATGCA 
    char Base[256];   // Base[256] = {'T','-','5','A','T','G','C','A'} (repeat 3 times);
    int diglen;       // diglen = 1 (strlen("5"))
    int count;        // count = 3
} PT;                 // Position Type

typedef struct {
    int dep;          // real depth count
    int pos;          // position of this mpileup line
    char ref;         // reference base
    int bc[5];        // base count
    char mk[5];       // marker
    char *base;       // read base
    int pt;           // postype count
    PT postype[100];
} RL;                 // Read Line

typedef struct {
    int pos;          // position
    int cov, alt_c;   // count of coverage and alteration
    char chr[20];     // chromosome
    char alt[255];    // alteration sequene
    char type[4];     // type: (SNP; Ins; Del)
    char ref;         // reference base
} RP;                 // report information

typedef struct {
    int chr_l;
    int seq_start;
    int seq_l;
    char chr[20];
} Chr;


typedef struct {
    int help;
    float cutoff;
    char in[PATH_MAX];
    char out[PATH_MAX];
    char ref[PATH_MAX];
} ARG;

void Usage ( void );
/* Usage of the Program [ called by main() ] */

void ParseArgs (ARG *, int, char **);

void Progress ( int , int );

void ChrFai ( Chr *, char * );
/* Get index file of fasta [ Called by main() ]*/

char *ChrBuild ( Chr *, char *);
/* Build array for original reference [ Called by main() ]*/

DL GetDigit ( const char * );
/* To get the digit length of the insertion and deletion with the postion 
 * [ called by GetIn() and WriteFa() ] */

void DealInde ( RL *, int, char * );
/* Process the insertion and deletion and store them [ called by GetInde()] */

int GetInde ( char *, RL * );
/* Get insertion and deletion shift width [ called by GetBase() ] */

int GetUnit ( RL *, char *, RP *, Chr *, int *, float );
/* Calculate the percent of base type and decide which to choose [ called by GetBase()] */

char *GetBase ( RL *, RP *, Chr *, int *, char *, int *, float );
/* Iterate the mpileup position line [ called by main() ] */

void LineSplit ( char *, RL * );
/* Split the mpileup line by "\t" with strtok [ called by main() ] */

void WriteFa ( char *, char *, Chr * );
/* Write out the final result [ Called by main() ] */

int GetReport ( RL *, RP *, Chr *, int *, int, int );
/* Get the report infomation [ Called by GetUnit() ] */

void WriteReport ( RP *, int );
/* Write out the report info about the substituted base [ Called by main() ] */

FILE *X_fopen ( char *filename );
/* Open file safely */
