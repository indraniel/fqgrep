/* I N C L U D E S ***********************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <zlib.h>  
#include "kseq.h"
#include "bm.h"
#include "bitap.h"

/* D E F I N E S *************************************************************/
#define VERSION "0.1"
#define PRG_NAME "fqgrep"
#define FASTQ_FILENAME_MAX_LENGTH 1024
#define MAX_PATTERN_LENGTH 1024

/* D A T A    S T R U C T U R E S ********************************************/
typedef struct {
    int count;
    int color;
    int mismatches;
    char output_fastq[FASTQ_FILENAME_MAX_LENGTH];
    char search_pattern[MAX_PATTERN_LENGTH];
} options;

/* 
   declare the type of file handler and the read() function  
   as described here:
   http://lh3lh3.users.sourceforge.net/parsefastq.shtml
*/
KSEQ_INIT(gzFile, gzread)  

/* P R O T O T Y P E S *******************************************************/
void help_message(void);
void version_info(void);
int  process_options(int argc, char *argv[], options *opts);
void search_input_fastq_file(FILE *out_fp,
                             const char *input_fastq, 
                             const options opts);
void display_read(FILE *out_fp, 
                  const kseq_t *seq, 
                  const char *substr_position_start, 
                  const options opts);
char* substring(const char *str, size_t start, size_t len);

/* G L O B A L S *************************************************************/


/* M A I N *******************************************************************/
int main(int argc, char *argv[]) {

    int opt_idx, index;
    FILE *out_fp;
    char input_fastq[FASTQ_FILENAME_MAX_LENGTH] = { '\0' };
    options opts = { 0, 0, 0, {'\0'}, {'\0'} };

    opt_idx = process_options(argc, argv, &opts);

    if (opt_idx >= argc) {
        fprintf(stderr, "%s : %s\n",
                        PRG_NAME, "[err] specify FASTQ files to process!");
        exit(1);
    }

    /* setup the appropriate output file pointer */
    if ( !strlen(opts.output_fastq) ) {
        out_fp = stdout;
    }
    else {
        if ( (out_fp = fopen(opts.output_fastq, "w")) == NULL ) {
            fprintf(stderr, "%s : [err] Could not open '%s' for writing.\n",
                            PRG_NAME, opts.output_fastq);
            exit(1);
        }
    }
    
    /* the remaining command line arguments are FASTQ(s) to process */
    while (opt_idx < argc) {
        strncpy(input_fastq, argv[opt_idx], FASTQ_FILENAME_MAX_LENGTH);
        search_input_fastq_file(out_fp, input_fastq, opts);
        opt_idx++;
    }

    fclose(out_fp);

    return 0;
}

/* F U N C T I O N S *********************************************************/
void 
help_message() {
    fprintf(stdout, "Usage: %s %s %s %s\n", 
                    PRG_NAME, "[options]", "-p <pattern>", "<fastq_files>");
    fprintf(stdout, "\t%-20s%-20s\n", "-h", "This help message");
    fprintf(stdout, "\t%-20s%-20s\n", "-v", "Program and version information");
    fprintf(stdout, "\t%-20s%-20s\n", "-p", "Pattern of interest to grep [REQUIRED]");
    fprintf(stdout, "\t%-20s%-20s\n", "-c", "Highlight matching string with color");
    fprintf(stdout, "\t%-20s%-20s\n", "-m", "Number of mismatches to allow for search pattern");
    fprintf(stdout, "\t%-20s%-20s\n", "", "[Default: 0]");
    fprintf(stdout, "\t%-20s%-20s\n", "-C", "Display only a total count of matches");
    fprintf(stdout, "\t%-20s%-20s\n", "", "(per input FASTQ file)");
    fprintf(stdout, "\t%-20s%-20s\n", "-o <fastq_file>", "Desired fastq output file.");
    fprintf(stdout, "\t%-20s%-20s\n", "", "If not specified, defaults to stdout");
}

void
version_info() {
    fprintf(stdout, "%s -- version: %s\n", PRG_NAME, VERSION);
}

int
process_options(int argc, char *argv[], options *opts) {
    int c;
    int index;
    char *opt_o_value = NULL;
    char *opt_p_value = NULL;

    while( (c = getopt(argc, argv, "hvcm:o:p:C")) != -1 ) {
        switch(c) {
            case 'h':
                help_message();
                exit(0);
                break;
            case 'v':
                version_info();
                exit(0);
                break;
            case 'o':
                opt_o_value = optarg;
                break;
            case 'p':
                opt_p_value = optarg;
                break;
            case 'c':
                opts->color = 1;
                break;
            case 'm':
                opts->mismatches = atoi(optarg);
                break;
            case 'C':
                opts->count = 1;
                break;
            case '?':
                exit(1);
             default:
                abort();
        }
    }

    /* ascertain whether a query pattern was given */
    if ( opt_p_value == NULL ) {
        fprintf(stderr, "%s : %s\n", PRG_NAME,
                        "[err] Specify a search pattern via the '-p' option!");
        fprintf(stderr, "Type '%s -h' for usage.\n", PRG_NAME);
        exit(1);
    }
    else {
        strncpy(opts->search_pattern, opt_p_value, MAX_PATTERN_LENGTH);
    }

    /* 
       ensure the number of mismatches does not exceed the search
       pattern length 
    */
    if ( opts->mismatches > strlen(opts->search_pattern) ) {
        fprintf(stderr, "%s : %s%d%s%zd%s\n", 
                        PRG_NAME,
                        "[err] The number of allowed mismatches (",
                        opts->mismatches,
                        ") exceeds the search pattern length (",
                        strlen(opts->search_pattern),
                        ")" );
        exit(1);
    }

    /* appropriately note output file (if given) */
    if ( opt_o_value != NULL ) {
        strncpy(opts->output_fastq, opt_o_value, FASTQ_FILENAME_MAX_LENGTH);
    }

    return optind;
}

void
search_input_fastq_file(FILE *out_fp, 
                        const char *input_fastq, 
                        const options opts) {
    gzFile fp;  
    kseq_t *seq;  
    int l, match_counter = 0;
    char *substr_position_start;

    // open the file handler  
    if ( (fp = gzopen(input_fastq, "r")) == NULL ) {
        fprintf(stderr, "%s : [err] Could not open FASTQ '%s' for reading.\n",
                        PRG_NAME, input_fastq);
        exit(1);
    } 

    // initialize seq  
    seq = kseq_init(fp); 

    // read sequence  
    while ( (l = kseq_read(seq)) >= 0 ) { 
        if ( strlen(seq->seq.s) < strlen(opts.search_pattern) ) {
            fprintf(stderr, "%s : %s '%s' %s (%zd) %s (%zd).\n",
                            PRG_NAME, 
                            "[err] For sequence ",
                            seq->name.s,
                            "search pattern length",
                            strlen(opts.search_pattern),
                            "exceeds sequence length",
                            strlen(seq->seq.s) );
            exit(1);
        }

        if (opts.mismatches == 0) {
            substr_position_start = 
                (char *) boyermoore_search( seq->seq.s, opts.search_pattern );
        }
        else {
            substr_position_start = 
                (char *) bitap_fuzzy_bitwise_search(
                        seq->seq.s, 
                        opts.search_pattern, 
                        opts.mismatches 
                );
        }

        if (substr_position_start != NULL) {
            match_counter++;
            if (opts.count == 0)
                display_read(out_fp, seq, substr_position_start, opts);
        }
    }

    kseq_destroy(seq); // destroy seq  
    gzclose(fp);       // close the file handler  

    //fprintf(stdout, "Mismatch param is %d\n", opts.mismatches);
    if (opts.count == 1) {
        if (match_counter == 1) {
            fprintf(out_fp, "%s : %d match\n", input_fastq, match_counter);
        }
        else {
            fprintf(out_fp, "%s : %d matches\n", input_fastq, match_counter);
        }
    }
}

void 
display_read(FILE *out_fp, 
             const kseq_t *seq, 
             const char *substr_position_start, 
             const options opts) {

    /* header portion of FASTQ read record */
    fprintf(out_fp, "@%s\n", seq->name.s);

    /* sequence portion of FASTQ read record */
    if (opts.color == 1) {
        char *sequence = seq->seq.s;
        size_t start, length;

        if (sequence == substr_position_start) {
            char *highlight = 
                substring( sequence, 0, strlen(opts.search_pattern) );

            start  = strlen(opts.search_pattern);
            length = strlen(sequence) - strlen(opts.search_pattern);

            char *remainder = 
                substring( sequence, start, length );

            fprintf(out_fp, "\033[31m%s\033[0m%s\n", highlight, remainder);
            free(highlight);
            free(remainder);
        }
        else {
            start  = 0;
            length = substr_position_start - sequence;
            char *begin = substring( sequence, start, length ); 

            start  = start + length;
            length = strlen(opts.search_pattern);
            char *highlight = substring( sequence, start, length );

            start  = start + length;
            length = strlen(sequence) - start;
            char *end = substring ( sequence, start, length );

            fprintf(out_fp, "%s\033[31m%s\033[0m%s\n", begin, highlight, end);
            free(begin);
            free(highlight);
            free(end);
        }
    }
    else {
        fprintf(out_fp, "%s\n", seq->seq.s);
    }

    /* comment portion of FASTQ read record */
    if (seq->comment.l) {
        fprintf(out_fp, "+%s\n", seq->comment.s);
    }
    else {
        fprintf(out_fp, "+\n");
    }

    /* quality portion of FASTQ read record */
    if (seq->qual.l) {
        fprintf(out_fp, "%s\n", seq->qual.s);
    }
    else {
        fprintf(out_fp, "\n", seq->qual.s);
    }
}

char* 
substring(const char *str, size_t start, size_t len) {
    if (str == NULL 
        || strlen(str) == 0 
        || strlen(str) < start 
        || strlen(str) < (start+len))
    return NULL;

  return strndup(str + start, len);
}
