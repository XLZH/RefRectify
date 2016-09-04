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

int main (int argc, char **argv)
{
    ARG Args; memset(&Args, 0, sizeof(ARG));
    ParseArgs(&Args, argc, argv);
    if ( Args.help || argc < 6 ) {
        Usage (); exit(1);
    }
    if ( Args.cutoff == 0 ) Args.cutoff = 0.9;

    char *outfa, *pstart, *pbase, *refer;
    RL readline; RP report[MEM_UNIT]; Chr chrom;
    int baselen, pre_pos, report_c;
    time_t start, end;

    baselen = pre_pos = report_c = 0;

    time (&start);
    ChrFai ( &chrom, Args.ref ); refer = ChrBuild ( &chrom, Args.ref );

    pstart = outfa = (char *)malloc((chrom.chr_l + 64) * sizeof(char) );
    FILE *mpileup = X_fopen ( Args.in );
    printf("[*] Start to go through the whole mpileup file......\n");

    while ( fgets(Buff, MAX_DEPTH, mpileup) )
    {
        LineSplit ( Buff, &readline );
        if (readline.pos % 2000 == 0 ) Progress ( readline.pos, chrom.chr_l );
        pbase = GetBase(&readline, report, &chrom, &report_c, refer, &pre_pos, Args.cutoff); 
        baselen = strlen(pbase); strcpy(outfa, pbase); outfa += baselen;
        free (readline.base); free(pbase);
    } *outfa = '\0'; printf("\n");
    WriteFa ( pstart, Args.out, &chrom ); WriteReport ( report, report_c );

    free (pstart); free(refer); fclose (mpileup); time(&end);
    printf ("Time Consum is: %.2lf (s)\n", difftime(end, start));
}


void Usage ( void )
{
    char *Usage =
    "\nProgram: RefRectify ( Rectify the reference [germline SNP, Insert, Delete])\n" \
    "Version: 3\n" \
    "\nUsage: RefRctify [options]\n\n" \
    "Options:\n" \
    "         [-h]        This will print helpInfo\n" \
    "         [-i]        Input file [mplieup file]\n" \
    "         [-o]        Output file [rectified reference]\n" \
    "         [-r]        Original reference that will be rectified\n"
    "         [-c]        Cutoff for alt-frequence(should bigger than 0.5) [default = 0.9]\n\n";
    fprintf( stderr, "%s", Usage );
}


void ParseArgs ( ARG* Args, int argc, char *argv[] )
{
    int ch;
    opterr = 0;

    while ((ch = getopt(argc, argv, "hi:o:r:c:")) != -1)
    {
        switch (ch) {
            case 'h': 
                Args->help = 1; break;
            case 'i': 
                strcpy(Args->in, optarg); break;
            case 'o': 
                strcpy(Args->out, optarg); break;
            case 'r': 
                strcpy(Args->ref, optarg); break;
            case 'c':
                Args->cutoff = atof(optarg); 
                if ( Args->cutoff <= 0.5 ) {
                    fprintf(stderr, "[Err::%s::%d] Cutoff should biger than 0.5!\n", __func__, __LINE__); exit(-1);
                } break;
            case '?':
                fprintf(stderr, "[Err::%s::%d] Option error occour!.\n", __func__, __LINE__); 
                Args->help = 1;
        }
    }
}


void Progress ( int pos, int all )
{
    float percent = (float)pos / all * 100;
    printf("\r");
    for (size_t i=0; i < 50; ++i ) {
        if ( i <= percent / 2 ) 
            printf("▆");
        else 
            printf("_");
    }
    printf(" (%.1f %%)",percent); fflush(stdout);
}


FILE *X_fopen ( char *filename )
{
    FILE *fp = fopen( filename, "r" );
    if ( !fp ) {
        fprintf( stderr, "[Err::%s::%d] Failed to open %s\n", __func__, __LINE__, filename ); exit (-1);
    }
    setvbuf (fp, NULL, _IOFBF, BUFSIZ);
    return fp;
}


void ChrFai ( Chr *chrom, char *fafile )
{
    char *pword; 
    char buff[512];
    int *start = &(chrom->chr_l);
    
    if ( !fafile ) {
        fprintf(stderr, "[Err::%s::%d] Can't find fasta index(.fai)!\n", __func__, __LINE__); exit(-1);
    }
    strcpy(buff, fafile); strcat (buff, ".fai");

    FILE *fp = X_fopen( buff );
    fgets(buff, 512, fp); pword = strtok (buff, "\t"); strcpy(chrom->chr, pword);
    for ( int i=0; i < 3 && pword; ++i ) {
        pword = strtok(NULL, "\t"); *start = atoi(pword); ++start;
    } fclose(fp);
}

char *ChrBuild ( Chr *chrom, char *fafile )
{
    FILE *fp; 
    char buff[256];
    char *fa, *current;

    if ( !(fa = current = (char *)malloc((chrom->chr_l + 64) * sizeof(char))) ) {
        fprintf(stderr, "[Err::%s::%d] Failed to allocat memory.\n", __func__, __LINE__); exit(-1);
    }
    fp = X_fopen( fafile ); fseek( fp, (long)chrom->seq_start, SEEK_SET );
    printf("[*] Start to build the original reference array......\n");
    while ( fgets(buff, 256, fp) ) {
        for (int i=0; i < chrom->seq_l; ++i ) {
            *current = toupper(buff[i]); ++current;
        }
    } *current = '\0';
    fclose(fp); return (fa);
}


DL GetDigit ( const char *inseq )
{
    char digit[5] = {'\0'};
    char *pdigit = digit;
    DL seqdig;

    for ( ++inseq; isdigit(*inseq); ++inseq )
        *pdigit++ = *inseq;
    seqdig.digit = atoi(digit); seqdig.diglen = strlen(digit);
    return seqdig;
}

int GetReport ( RL *readline, RP *report, Chr *chrom, int *report_c, int type, int index )
{
    int curr_index = *report_c;

    strcpy(report[curr_index].chr, chrom->chr);
    report[curr_index].pos = readline->pos; report[curr_index].ref = readline->ref;
    report[curr_index].cov = readline->dep;
    if ( type == 1 ) {
        report[curr_index].alt_c = readline->bc[index];
        report[curr_index].alt[0] = readline->mk[index]; report[curr_index].alt[1] = '\0';
        strcpy(report[curr_index].type, "SNP"); (*report_c)++; return 0;
    } 

    report[curr_index].alt_c = readline->postype[index].count;
    strcpy(report[curr_index].alt, readline->postype[index].Base);
    if ( type == 2 ) {
        strcpy(report[curr_index].type, "INS"); (*report_c)++; return 0;
    }
    if ( type == 3 ) {
       strcpy(report[curr_index].type, "DEL"); (*report_c)++; return 0;
    }
}

int GetUnit ( RL *readline, char *outunit, RP *report, Chr *chrom, int *report_c, float cutoff )
{
    if ( readline->dep >=0 && readline->dep < 10) {
        *outunit = readline->ref; return 1;
    }
    if ( (float)readline->bc[0] / readline->dep > cutoff ) {
        *outunit = readline->ref; return 0;
    }
    for ( int i=1; i < 5; ++i ) {
        if ( (float)readline->bc[i] / readline->dep > cutoff ) {
            *outunit = readline->mk[i]; 
            GetReport(readline, report, chrom, report_c, 1, i); return 0;
        }
    }
    for ( int i=0; i < readline->pt; ++i ) {
        if ( (float)readline->postype[i].count / readline->dep > cutoff ) {
            if ( readline->postype[i].Base[1] == 43 ) {
                strcpy (outunit, readline->postype[i].Base); 
                GetReport (readline, report, chrom, report_c, 2, i); return 0;
            } else {
                strncpy (outunit, readline->postype[i].Base, (readline->postype[i].diglen + 2));
                GetReport (readline, report, chrom, report_c, 3, i); return 0;
            }
        }
    }
    *outunit = readline->ref; return 1;
}


void DealInde ( RL *readline, int digit, char *inde)
{
    int exist = 0;
    for (int i=0; i < readline->pt; ++i) {
        if ( !(strcmp(inde, readline->postype[i].Base)) ) {
            ++readline->postype[i].count; exist = 1; break;
        }
    }
    if ( !exist ) {
        strcpy(readline->postype[readline->pt].Base, inde);
        ++readline->postype[readline->pt].count;
        readline->postype[readline->pt].diglen = digit; 
        ++readline->pt;
    }
    if ( inde[0] == readline->ref ) --readline->bc[0];
    else {
        for ( int i=1; i < 5; ++i ) {
            if ( inde[i] == readline->mk[i] ) { --readline->bc[i]; break; }
        }
    }
}


int GetInde ( char *inseq, RL *readline )
{
    char inde[256] = {'\0'}; 
    char *pinde = &inde[1];
    DL seqdig = GetDigit ( inseq );
    int indelen = seqdig.diglen + seqdig.digit + 1;

    inde[0] = TRANS(*(inseq-1), readline->ref);

    for ( int i=1; i <= indelen; ++i,++pinde,++inseq) 
    {
        if (i == 1) *pinde = *inseq;
        else *pinde = toupper(*inseq);
    }
    DealInde (readline, seqdig.diglen, inde);
    return (indelen - 1);
}


char *GetBase ( RL *readline, RP *report, Chr *chrom, int *report_c, char *fa, int *pre_pos, float cutoff) // Core part
{
    char *b = readline->base; 
    char *start, *outunit;
    int misspos_c;

    if ( readline->pos - *pre_pos != 1 ) {
        misspos_c = readline->pos - *pre_pos - 1;
        start = outunit = (char *) calloc( (misspos_c + 256), sizeof(char) );
        strncpy( outunit, (fa + *pre_pos), misspos_c );
        outunit += (size_t)misspos_c; *pre_pos = readline->pos;
    } 
    else {
        start = outunit = (char *) calloc( 256, sizeof(char) );
        *pre_pos = readline->pos;
    }
    for ( b; *b != '\0'; ++b )
    {
        if ( *b == 46 || *b == 44 ) {  // [ . , ]
            ++readline->dep; ++readline->bc[0]; continue; 
        }
        if ( *b == 36 || *b == 42 || *b == 94 || *b == 93) continue; // [ $ * ] ^ ]
        for ( int i=1; i < 5; ++i) { 
            if (toupper(*b) == readline->mk[i]) {
                ++readline->dep; ++readline->bc[i]; break;
            }
        }
        if ( *b == 43 || *b == 45 ) b += GetInde (b, readline); // [ + - ]
    }
    GetUnit (readline, outunit, report, chrom, report_c, cutoff); 
    return start;
}


void LineSplit ( char* bufline, RL *readline )
{
    memset(readline, 0, sizeof(RL));
    readline->base = (char *) malloc(strlen(bufline) * sizeof(char));
    char *pword = strtok(bufline, "\t");

    for ( int i=0; pword != NULL; ++i )
    {
        pword = strtok(NULL, "\t");
        if ( i == 0 ) { readline->pos = atoi(pword); continue; }
        if ( i == 1 ) { readline->ref = toupper(*pword); continue; }
        if ( i == 3 ) { strcpy(readline->base, pword); break; }
    }
    strcpy (readline->mk, "!ATGC"); // "!" is a dummy base.
}


void WriteFa ( char *start, char *outfile, Chr *chrom)
{
    FILE *Fasta = fopen(outfile, "a");
    DL marklen;
    
    printf("[*] Start to output the rectified reference......\n");
    fprintf(Fasta, ">%s\t%s\n", chrom->chr, "Complete Genome");
    while ( *start != '\0' ) {
        for ( int i=0; i < 50 && *start != '\0'; ++i ) {
            if ( *start == 0x2b ) { // [ + ]
                marklen = GetDigit ( start ); 
                start += (marklen.diglen + 1); --i;
            } 
            else if ( *start == 0x2d ) { // [ - ]
                marklen = GetDigit ( start ); 
                start += (marklen.diglen + marklen.digit + 1); --i;
            }
            else { // [ A T G C ]
                fputc (*start, Fasta); ++start;
            }
        } fputc ('\n', Fasta);
    }
    fclose (Fasta);
}

void WriteReport ( RP *report, int report_c )
{
    FILE *reportfile = fopen ("RefRectify.report","w");
    
    printf("[*] Start to output the infomation of rectified position......\n");
    fprintf(reportfile, "CHROM\tPOS\tTYPE\tREF\tALT\tCOVERAGE:ALTCOUNT\n");
    for ( int i=0; i < report_c; ++i ) {
        fprintf(reportfile, "%s\t%d\t%s\t%c\t%s\t%d:%d\n", report[i].chr, report[i].pos, report[i].type, report[i].ref, report[i].alt, report[i].cov, report[i].alt_c);
    }
    fclose(reportfile);
}
