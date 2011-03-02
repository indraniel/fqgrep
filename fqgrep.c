/* L I C E N S E *************************************************************/

/*
    Copyright (C) 2010, 2011 Indraniel Das <indraniel@gmail.com>
                             and Washington University in St. Louis

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, see <http://www.gnu.org/licenses/>
*/

/* I N C L U D E S ***********************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <string.h>
#include <unistd.h>
#include <zlib.h>  
#include <tre/tre.h>
#include "kseq.h"
#include "bm.h"

/* D E F I N E S *************************************************************/
#define VERSION "0.3"
#define PRG_NAME "fqgrep"
#define FASTQ_FILENAME_MAX_LENGTH 1024
#define MAX_PATTERN_LENGTH 1024
#define MAX_DELIM_LENGTH 10

/* D A T A    S T R U C T U R E S ********************************************/
typedef struct {
    int count;
    int color;
    int force_tre;
    int report_fastq;
    int report_fasta;
    int report_stats;
    int max_mismatches;
    int cost_insertions;
    int cost_deletions;
    int cost_substitutions;
    int max_insertions;
    int max_deletions;
    int max_substitutions;
    char output_fastq[FASTQ_FILENAME_MAX_LENGTH];
    char search_pattern[MAX_PATTERN_LENGTH];
    char delim[MAX_DELIM_LENGTH];         /* delimiter used in stats report */
    regex_t *tre_regex;                   /* Compiled tre regexp */
    regaparams_t *tre_regex_match_params; /* tre regexp matching parameters */
} options;

typedef struct {
    char *sequence;
    char *substr_start;
    char *substr_end;
    int  start_pos;
    int  end_pos;
    int  num_mismatches;
    int  num_insertions;
    int  num_deletions;
    int  num_substitutions;
} read_match;

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
void  report_read(FILE *out_fp,
                  const options *opts,
                  const kseq_t *seq,
                  const read_match *info);
void  report_fastq(FILE *out_fp,
                   const options *opts,
                   const kseq_t *seq,
                   const read_match *info);
void  report_fasta(FILE *out_fp,
                   const options *opts,
                   const kseq_t *seq,
                   const read_match *info);
void  report_stats(FILE *out_fp,
                   const options *opts,
                   const kseq_t *seq,
                   const read_match *info);
void  display_sequence(FILE *out_fp,
                       const options *opts,
                       const char *sequence,
                       const char *substr_start,
                       const char *substr_end,
                       const int  start_pos,
                       const int  end_pos);
void  setup_tre(regaparams_t *params, regex_t *regexp, options *opts);
void  approximate_regexp_search(const options *opts, read_match *info);
char* substring(const char *str, size_t start, size_t len);
char* stringn_duplicate(const char *str, size_t n);

/* G L O B A L S *************************************************************/


/* M A I N *******************************************************************/
int main(int argc, char *argv[]) {

    int opt_idx, index;
    FILE *out_fp;
    char input_fastq[FASTQ_FILENAME_MAX_LENGTH] = { '\0' };

    /* application of default options */
    options opts = {
        0,            // count flag
        0,            // color flag
        0,            // force tre engine flag
        1,            // output fastq report
        0,            // output fasta report
        0,            // output stats report
        0,            // max mismatches allowed
        1,            // cost of insertions
        1,            // cost of deletions
        1,            // cost of substitutions
        INT_MAX,      // maxiumum allowable insertions in match
        INT_MAX,      // maxiumum allowable deletions in match
        INT_MAX,      // maxiumum allowable substitutions in match
        {'\0'},       // output fastq file name
        {'\0'},       // search pattern string
        "\t",         // delimiter string for stats report
        NULL,         // pointer to tre regexp entity
        NULL          // pointer to tre regexp matching parameters
    };

    opt_idx = process_options(argc, argv, &opts);

    if (opt_idx >= argc) {
        fprintf(stderr, "%s : %s\n",
                        PRG_NAME, "[err] specify FASTQ files to process!");
        exit(1);
    }

    /* setup and compile the tre regexp if needed */
    if (opts.max_mismatches != 0 || opts.force_tre == 1) {
        regex_t regxp;                    /* Compiled pattern to search for. */
        regaparams_t match_params;        /* regexp matching parameters */

        setup_tre( &match_params, &regxp, &opts );

//    fprintf(stdout, "TRE regex params setup:\n");
//    fprintf(stdout, "\t%-12s : %4d\n", "cost_ins",   match_params.cost_ins);
//    fprintf(stdout, "\t%-12s : %4d\n", "cost_del",   match_params.cost_del);
//    fprintf(stdout, "\t%-12s : %4d\n", "cost_subst", match_params.cost_subst);
//    fprintf(stdout, "\t%-12s : %4d\n", "max_cost",   match_params.max_cost);

//    fprintf(stdout, "\t%-12s : %4d\n", "max_ins",   match_params.max_ins);
//    fprintf(stdout, "\t%-12s : %4d\n", "max_del",   match_params.max_del);
//    fprintf(stdout, "\t%-12s : %4d\n", "max_subst", match_params.max_subst);
//    fprintf(stdout, "\t%-12s : %4d\n", "max_err",   match_params.max_err);
//    fprintf(stdout, "\n\n");
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
    fprintf(stdout, "\t%-20s%-20s\n", "-p <STRING>", "Pattern of interest to grep [REQUIRED]");
    fprintf(stdout, "\t%-20s%-20s\n", "-c", "Highlight matching string with color");
    fprintf(stdout, "\t%-20s%-20s\n", "-f", "Output matches in FASTA format");
    fprintf(stdout, "\t%-20s%-20s\n", "-r", "Output matches in detailed stats report format");
    fprintf(stdout, "\t%-20s%-20s\n", "-b <STRING>", "Delimiter string to separate columns");
    fprintf(stdout, "\t%-20s%-20s\n", "", "in detailed stats report [Default: '\\t']");
    fprintf(stdout, "\t%-20s%-20s\n", "-m <INT>", "Total number of mismatches to at most allow for");
    fprintf(stdout, "\t%-20s%-20s\n", "", "in search pattern [Default: 0]");
    fprintf(stdout, "\t%-20s%-20s\n", "-s <INT>", "Max threshold of substitution mismatches to allow");
    fprintf(stdout, "\t%-20s%-20s\n", "", "for in search pattern [Default: unlimited]");
    fprintf(stdout, "\t%-20s%-20s\n", "-i <INT>", "Max threshold of insertion mismatches to allow for");
    fprintf(stdout, "\t%-20s%-20s\n", "", "in search pattern [Default: unlimited]");
    fprintf(stdout, "\t%-20s%-20s\n", "-d <INT>", "Max threshold of deletion mismatches to allow for");
    fprintf(stdout, "\t%-20s%-20s\n", "", "in search pattern [Default: unlimited]");
    fprintf(stdout, "\t%-20s%-20s\n", "-S <INT>", "Cost of base substitutions in obtaining");
    fprintf(stdout, "\t%-20s%-20s\n", "", "approximate match [Default: 1]");
    fprintf(stdout, "\t%-20s%-20s\n", "-I <INT>", "Cost of base insertions in obtaining");
    fprintf(stdout, "\t%-20s%-20s\n", "", "approximate match [Default: 1]");
    fprintf(stdout, "\t%-20s%-20s\n", "-D <INT>", "Cost of base deletions in obtaining");
    fprintf(stdout, "\t%-20s%-20s\n", "", "approximate match [Default: 1]");
    fprintf(stdout, "\t%-20s%-20s\n", "-e", "Force tre regexp engine usage");
    fprintf(stdout, "\t%-20s%-20s\n", "-C", "Display only a total count of matches");
    fprintf(stdout, "\t%-20s%-20s\n", "", "(per input FASTQ file)");
    fprintf(stdout, "\t%-20s%-20s\n", "-o <out_file>", "Desired output file.");
    fprintf(stdout, "\t%-20s%-20s\n", "", "If not specified, defaults to stdout");
}

void
version_info() {
    fprintf(stdout, "%s -- version: %s\n", PRG_NAME, VERSION);
    fprintf(stdout, "%s --  %s:%s\n",
            "Author", "Indraniel Das",
            "<idas at wustl dot edu> or <indraniel at gmail dot com>");
}

int
process_options(int argc, char *argv[], options *opts) {
    int c;
    int index;
    char *opt_o_value = NULL;
    char *opt_p_value = NULL;
    char *opt_b_value = NULL;

    while( (c = getopt(argc, argv, "hvecfrm:i:s:d:o:p:b:CD:I:S:")) != -1 ) {
        switch(c) {
            case 'h':
                help_message();
                exit(0);
                break;
            case 'v':
                version_info();
                exit(0);
                break;
            case 'e':
                opts->force_tre = 1;
                break;
            case 'o':
                opt_o_value = optarg;
                break;
            case 'p':
                opt_p_value = optarg;
                break;
            case 'b':
                opt_b_value = optarg;
                break;
            case 'c':
                opts->color = 1;
                break;
            case 'f':
                opts->report_fasta = 1;
                opts->report_fastq = 0;
                opts->report_stats = 0;
                break;
            case 'r':
                opts->report_fasta = 0;
                opts->report_fastq = 0;
                opts->report_stats = 1;
                break;
            case 'm':
                opts->max_mismatches = atoi(optarg);
                break;
            case 'i':
                opts->max_insertions = atoi(optarg);
                break;
            case 's':
                opts->max_substitutions = atoi(optarg);
                break;
            case 'd':
                opts->max_deletions = atoi(optarg);
                break;
            case 'I':
                opts->cost_insertions = atoi(optarg);
                break;
            case 'D':
                opts->cost_deletions = atoi(optarg);
                break;
            case 'S':
                opts->cost_substitutions = atoi(optarg);
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

    /* setup delimiter for stats report (if given) */
    if ( opt_b_value != NULL ) {
        strncpy(opts->delim, opt_b_value, MAX_DELIM_LENGTH);
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
    read_match match_info;

    // open the file handler
    if ( strcmp(input_fastq, "-") == 0 ) {
        fp = gzdopen(fileno(stdin), "r");
    }
    else {
        fp = gzopen(input_fastq, "r");
    }

    if ( (fp == NULL) && (strcmp(input_fastq, "-") != 0) ) {
        fprintf(stderr, "%s : [err] Could not open FASTQ '%s' for reading.\n",
                        PRG_NAME, input_fastq);
        exit(1);
    }

    if ( (fp == NULL) && (strcmp(input_fastq, "-") == 0) ) {
        fprintf(stderr, "%s : [err] Could not open stdin for reading.\n",
                        PRG_NAME);
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

        /* initialize the match info structure */
        match_info.sequence     = seq->seq.s;
        match_info.substr_start = NULL;
        match_info.substr_end   = NULL;
        match_info.start_pos    = 0;
        match_info.end_pos      = 0;

        match_info.num_mismatches    = 0;
        match_info.num_insertions    = 0;
        match_info.num_deletions     = 0;
        match_info.num_substitutions = 0;

        if (opts.max_mismatches == 0 && opts.force_tre == 0) {
//            fprintf(stdout, "Running boyer moore search\n");
            match_info.substr_start =
                (char *) boyermoore_search( seq->seq.s, opts.search_pattern );
            if (match_info.substr_start != NULL) {
                match_info.substr_end =
                    match_info.substr_start + strlen(opts.search_pattern);
                match_info.start_pos =
                    (int) (match_info.substr_start - seq->seq.s);
                match_info.end_pos =
                    (int) ( match_info.start_pos + strlen(opts.search_pattern) );
            }
        }
        else {
//            fprintf(stdout, "Running TRE search\n");
            approximate_regexp_search( &opts, &match_info );
        }

        if (match_info.substr_start != NULL) {
            match_counter++;
            if (opts.count == 0)
                report_read( out_fp, &opts, seq, &match_info );
        }
    }

    kseq_destroy(seq); // destroy seq  
    gzclose(fp);       // close the file handler  

    //fprintf(stdout, "Mismatch param is %d\n", opts.max_mismatches);
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
report_read(FILE *out_fp,
            const options *opts,
            const kseq_t *seq,
            const read_match *info) {
    if (opts->report_fasta) {
        report_fasta(out_fp, opts, seq, info);
    }
    else if (opts->report_stats) {
        report_stats(out_fp, opts, seq, info);
    }
    else {
        report_fastq(out_fp, opts, seq, info);
    }
}

void
display_sequence(FILE *out_fp,
                 const options *opts,
                 const char *sequence,
                 const char *substr_start,
                 const char *substr_end,
                 const int  start_pos,
                 const int  end_pos) {
    if (opts->color == 1) {
        size_t start, length;

        if (sequence == substr_start) {
            char *highlight = 
                substring( sequence, 0, substr_end - sequence );

            start  = substr_end - substr_start;
            length = strlen(sequence) -
                     ( substr_end - substr_start );

            char *remainder =
                substring( sequence, start, length );

            fprintf(out_fp, "\033[31m%s\033[0m%s", highlight, remainder);
            free(highlight);
            free(remainder);
        }
        else {
            start  = 0;
            length = substr_start - sequence;
            char *begin = substring( sequence, start, length );

            start  = start + length;
            length = (size_t) (end_pos - start_pos);
            char *highlight = substring( sequence, start, length );

            start  = start + length;
            length = strlen(sequence) - start;
            char *end = substring( sequence, start, length );

            fprintf(out_fp, "%s\033[31m%s\033[0m%s", begin, highlight, end);
            free(begin);
            free(highlight);
            free(end);
        }
    }
    else {
        fprintf(out_fp, "%s", sequence);
    }
}

void
report_stats(FILE *out_fp,
             const options *opts,
             const kseq_t *seq,
             const read_match *info) {

    static int header_flag = 0;
    char *match = NULL;

    /*
       stat report columns are
         1. read_name
         2. mismatches
         3. num_ins
         4. num_del
         5. num_subst
         6. start_position
         7. end_position
         8. match string
         9. sequence string
        10. quality string (if available)
     */

    if (header_flag == 0) {
        fprintf(out_fp, "%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s",
                "read name",
                opts->delim,
                "total mismatches",
                opts->delim,
                "# insertions",
                opts->delim,
                "# deletions",
                opts->delim,
                "# substitutions",
                opts->delim,
                "start position",
                opts->delim,
                "end position",
                opts->delim,
                "match string",
                opts->delim,
                "sequence"
        );

         /* quality string portion of header */
        if (seq->qual.l) {
            fprintf(out_fp, "%s", opts->delim);
            fprintf(out_fp, "%s", "quality");
        }
        fprintf(out_fp, "\n");
        header_flag = 1;
    }

    fprintf(out_fp, "%s%s%d%s%d%s%d%s%d%s%d%s%d%s",
            seq->name.s,
            opts->delim,
            info->num_mismatches,
            opts->delim,
            info->num_insertions,
            opts->delim,
            info->num_deletions,
            opts->delim,
            info->num_substitutions,
            opts->delim,
            info->start_pos,
            opts->delim,
            info->end_pos,
            opts->delim
    );

    /* match string portion of stats report */
    size_t start, length;

    start  = (size_t) info->start_pos;
    length = (size_t) (info->end_pos - info->start_pos);

    match = substring( seq->seq.s, start, length );
    fprintf(out_fp, "%s%s", match, opts->delim);
    free(match);

    /* sequence portion of stats report */
    display_sequence (out_fp,
                      opts,
                      seq->seq.s,
                      info->substr_start,
                      info->substr_end,
                      info->start_pos,
                      info->end_pos);

    /* quality string portion of stats report */
    if (seq->qual.l) {
        fprintf(out_fp, "%s", opts->delim);
        fprintf(out_fp, "%s", seq->qual.s);
    }

    /* termination of record line */
    fprintf(out_fp, "\n");
}

void
report_fasta(FILE *out_fp,
             const options *opts,
             const kseq_t *seq,
             const read_match *info) {
    /* header portion of FASTA read record */
    fprintf(out_fp, ">%s\n", seq->name.s);

    /* sequence portion of FASTA read record */
    display_sequence (out_fp,
                      opts,
                      seq->seq.s,
                      info->substr_start,
                      info->substr_end,
                      info->start_pos,
                      info->end_pos);
    fprintf(out_fp, "\n");
}

void
report_fastq(FILE *out_fp,
             const options *opts,
             const kseq_t *seq,
             const read_match *info) {

    /* header portion of FASTQ read record */
    fprintf(out_fp, "@%s\n", seq->name.s);

    /* sequence portion of FASTQ read record */
    display_sequence (out_fp,
                      opts,
                      seq->seq.s,
                      info->substr_start,
                      info->substr_end,
                      info->start_pos,
                      info->end_pos);
    fprintf(out_fp, "\n");

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

void
setup_tre(regaparams_t *params, regex_t *regexp, options *opts) {
    /* Step 1: setup the TRE regexp matching parameters */

    /* setup the default match parameters */
    tre_regaparams_default(params);

    /* Set the maximum number of errors allowed for a record to match. */
    params->max_cost = opts->max_mismatches;

    /* Set the insertion, deletion, substitution costs. */
    params->cost_ins = opts->cost_insertions;
    params->cost_del = opts->cost_deletions;
    params->cost_subst = opts->cost_substitutions;

    /* Set the max insertion, deletion, substitution allowances. */
    params->max_ins = opts->max_insertions;
    params->max_del = opts->max_deletions;
    params->max_subst = opts->max_substitutions;

    /* Step 2: compile the regex */

    /*
       always allowing for POSIX extended regular expression syntax
       (REG_EXTENDED)
       always being case insenstive with the regular expression
       (REG_ICASE)
    */
    int errcode;
    int comp_flags  = REG_EXTENDED | REG_ICASE ;
    errcode = tre_regcomp(regexp, opts->search_pattern, comp_flags);
    if (errcode) {
        char errbuf[256];
        tre_regerror(errcode, regexp, errbuf, sizeof(errbuf));
        fprintf(stderr, "%s: %s: -- %s -- %s\n",
          PRG_NAME,
              "Error in compiling search pattern",
              opts->search_pattern,
              errbuf
        );
        exit(1);
    }

    /* Step 3: assign to opts hash */
    opts->tre_regex = regexp;
    opts->tre_regex_match_params = params;
}

void
approximate_regexp_search(const options *opts, read_match *info) {
    int errcode;
    regmatch_t pmatch = { 0, 0 };     /* matched pattern structure */
    regamatch_t match;                /* overall match structure */

    /* initialize the overall match struct */
    memset(&match, 0, sizeof(match));
    /* assign default pattern structure */
    match.pmatch = &pmatch;
    /* initialization of pmatch array */
    match.nmatch = 1;

    /* perform the regexp search on the sequence string */
    errcode = tre_regaexec(
            opts->tre_regex,
            info->sequence,
            &match,
            *(opts->tre_regex_match_params),
            0
    );

    if (errcode != REG_OK) {
//        fprintf(stdout, "Found no matches!\n");
//        fprintf(stdout, "%6s : %s\n", "regexp", opts->search_pattern);
//        fprintf(stdout, "%6s : %s\n", "record", sequence);
        return;
    }

    /* found a match! */

    info->num_mismatches    = match.cost;
    info->num_insertions    = match.num_ins;
    info->num_deletions     = match.num_del;
    info->num_substitutions = match.num_subst;
    info->start_pos         = pmatch.rm_so;
    info->end_pos           = pmatch.rm_eo;

    info->substr_start      = info->sequence + (size_t) info->start_pos;
    info->substr_end        = info->substr_start + (size_t) info->end_pos;

//    fprintf(stdout, "Found match!\n");
//    fprintf(stdout, "\t%10s : %s\n", "record", info->sequence);
//    fprintf(stdout, "\t%10s : %s\n", "pattern", opts->search_pattern);
//    fprintf(stdout, "\t%10s : %4d\n", "cost", match.cost);
//    fprintf(stdout, "\t%10s : %4d\n", "num_ins", match.num_ins);
//    fprintf(stdout, "\t%10s : %4d\n", "num_del", match.num_del);
//    fprintf(stdout, "\t%10s : %4d\n", "num_subst", match.num_subst);
//    fprintf(stdout, "\t%10s : %2d - %2d\n", "char pos",
//                                            pmatch.rm_so, pmatch.rm_eo);

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
