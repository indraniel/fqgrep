/* I N C L U D E S ***********************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <zlib.h>  
#include "kseq.h"
#include "bm.h"
#include "libbitap.h"

/* D E F I N E S *************************************************************/
#define VERSION "0.2"
#define PRG_NAME "fqgrep"
#define FASTQ_FILENAME_MAX_LENGTH 1024
#define MAX_PATTERN_LENGTH 1024

/* D A T A    S T R U C T U R E S ********************************************/
typedef struct {
    int count;
    int color;
    int force_bitap;
    int mismatches;
    char output_fastq[FASTQ_FILENAME_MAX_LENGTH];
    char search_pattern[MAX_PATTERN_LENGTH];
    bitapType *bitap_regexp;  // compiled libbitap regexp
} options;

/* 
   declare the type of file handler and the read() function
   as described here:
   http://lh3lh3.users.sourceforge.net/parsefastq.shtml
*/
KSEQ_INIT(gzFile, gzread)  

/* P R O T O T Y P E S *******************************************************/
void  help_message(void);
void  version_info(void);
int   process_options(int argc, char *argv[], options *opts);
void  search_input_fastq_file(FILE *out_fp,
                              const char *input_fastq,
                              const options opts);
void  display_read(FILE *out_fp,
                   const kseq_t *seq,
                   const char *substr_position_start,
                   const char *substr_position_end,
                   const options opts);
void  compile_bitap_regexp(bitapType *b, options *opts);
char* substring(const char *str, size_t start, size_t len);
char* stringn_duplicate(const char *str, size_t n);

/* G L O B A L S *************************************************************/


/* M A I N *******************************************************************/
int main(int argc, char *argv[]) {

    int opt_idx, index;
    FILE *out_fp;
    char input_fastq[FASTQ_FILENAME_MAX_LENGTH] = { '\0' };
    options opts = { 0, 0, 0, 0, {'\0'}, {'\0'}, NULL };

    opt_idx = process_options(argc, argv, &opts);

    if (opt_idx >= argc) {
        fprintf(stderr, "%s : %s\n",
                        PRG_NAME, "[err] specify FASTQ files to process!");
        exit(1);
    }

    /* compile bitap regexp if needed */
    if (opts.mismatches != 0 || opts.force_bitap == 1) {
        bitapType bitap;
        compile_bitap_regexp( &bitap, &opts );
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
    fprintf(stdout, "\t%-20s%-20s\n", "-b", "Force bitap algorithm usage");
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

    while( (c = getopt(argc, argv, "hvbcm:o:p:C")) != -1 ) {
        switch(c) {
            case 'h':
                help_message();
                exit(0);
                break;
            case 'v':
                version_info();
                exit(0);
                break;
            case 'b':
                opts->force_bitap = 1;
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
    char *substr_position_start = NULL;
    char *substr_position_end = NULL;
    int actual_mismatches = 0;

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

        if (opts.mismatches == 0 && opts.force_bitap == 0) {
//            fprintf(stdout, "Running boyer moore search\n");
            substr_position_start = 
                (char *) boyermoore_search( seq->seq.s, opts.search_pattern );
            if (substr_position_start != NULL) {
                substr_position_end =
                    substr_position_start + strlen(opts.search_pattern);
            }
        }
        else {
//            fprintf(stdout, "Running libbitap search\n");
            substr_position_end =
                 FindWithBitap(
                        opts.bitap_regexp,
                        seq->seq.s,
                        strlen(seq->seq.s),
                        opts.mismatches,
                        &actual_mismatches,
                        &substr_position_start
                );
//            fprintf(stdout, "substr_position_end : %p (%s) \n", substr_position_end, substr_position_end);
//            fprintf(stdout, "substr_position_start : %p (%s) \n", substr_position_start, substr_position_start);
//            fprintf(stdout, "mismatches : %d\n", opts.mismatches);
//            fprintf(stdout, "acutal mismatches : %d\n", actual_mismatches);
        }

        if (substr_position_start != NULL) {
            match_counter++;
            if (opts.count == 0)
                display_read(
                        out_fp,
                        seq,
                        substr_position_start,
                        substr_position_end,
                        opts
                );
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
             const char *substr_position_end,
             const options opts) {

    /* header portion of FASTQ read record */
    fprintf(out_fp, "@%s\n", seq->name.s);

    /* sequence portion of FASTQ read record */
    if (opts.color == 1) {
        char *sequence = seq->seq.s;
        size_t start, length;

        if (sequence == substr_position_start) {
            char *highlight = 
                substring( sequence, 0, substr_position_end - sequence );

            start  = substr_position_end - substr_position_start;
            length = strlen(sequence) -
                     ( substr_position_end - substr_position_start );

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
            length = substr_position_end - substr_position_start;
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

void compile_bitap_regexp(bitapType *b, options *opts) {
    size_t search_pattern_length;
    char *bitap_search_pattern;

    // We are implicitly looking for substring matches.  We do not
    // want to match against the whole input read. According to
    // http://rational.co.za/libbitap/ we must prefix and postfix the
    // search pattern with '.*' for this to occur properly.

    search_pattern_length = strlen(".*")
                            + strlen(opts->search_pattern)
                            + strlen(".*");

    // account for NULL termination
    bitap_search_pattern = malloc(search_pattern_length + 1);
    if ( bitap_search_pattern == NULL ) {
        fprintf(stderr, "%s : %s %s\n",
                        PRG_NAME,
                        "Trouble mallocing for bitap search pattern string.",
                        "Out of memory!");
        exit(1);
    }

    // initialize the string to NULLs
    memset(bitap_search_pattern, '\0', search_pattern_length + 1);

    // construct the appropriate bitap search pattern string
    strcat(bitap_search_pattern, ".*");
    strcat(bitap_search_pattern, opts->search_pattern);
    strcat(bitap_search_pattern, ".*");

//    fprintf(stdout, "bitap_search_pattern is: |%s|\n", bitap_search_pattern);

    // compile the bitap search string according to the libbitap protocol
    int result = NewBitap(b, bitap_search_pattern);
    if (result < 0) {
        fprintf(stderr, "%s : %s '%s' %s %d\n",
                        PRG_NAME,
                        "[err] Trouble compiling bitap search pattern",
                        opts->search_pattern,
                        "at character position ",
                        (-1) * result );
    }

    // assign it to our options data structure
    opts->bitap_regexp = b;

    // free the search string
    free(bitap_search_pattern);
}

char* 
substring(const char *str, size_t start, size_t len) {
    char *substr;

    if (str == NULL 
        || strlen(str) == 0 
        || strlen(str) < start 
        || strlen(str) < (start+len))
    return NULL;

    substr = stringn_duplicate(str + start, len);
    return substr;
}

/*
   'stringn_duplicate' is really a poor man's duplication of glibc's
   'strndup'. However not all types of UNIXes implement strndup (like
   Apple, & FreeBSD). So the function below is used instead for the ease
   of portability.

   The implementation is taken from :
     http://webmail.solomo.de/~flo/strndup.patch
     http://unix.derkeiler.com/Mailing-Lists/FreeBSD/arch/2008-12/msg00010.html
*/
char*
stringn_duplicate(const char *str, size_t n) {
    size_t len;
    char *copy;

    for (len = 0; len < n && str[len]; len++)
        continue;

    if ( ( copy = malloc(len + 1) ) == NULL ) {
        fprintf(stderr, "%s : %s\n",
                        PRG_NAME, "Trouble with malloc. Out of memory!");
        exit(1);
    }

    memcpy(copy, str, len);
    copy[len] = '\0';
    return copy;
}
