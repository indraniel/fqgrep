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

/* N O T E S *****************************************************************/
/* compliation line:
   gcc -g -o tre-test -ltre tre-test.c

   This is just a simple program to test out and understand
   how to use the libtre library (http://laurikari.net/tre/).
*/

/* I N C L U D E S ***********************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <string.h>
#include <unistd.h>
#include <tre/tre.h>

/* D E F I N E S *************************************************************/
#define PRG_NAME "tre-test"
#define MAX_CHAR 1024

/* P R O T O T Y P E S *******************************************************/
int  process_options(int argc, char *argv[]);
void help_message(void);

/* G L O B A L S *************************************************************/
int cost_ins  = 1;
int cost_del  = 1;
int cost_subs = 1;
int max_cost  = 1;
char record[MAX_CHAR] = "";
char regexp[MAX_CHAR] = "";

/* M A I N *******************************************************************/
int main(int argc, char *argv[]) {
    int opt_index;

    opt_index = process_options(argc, argv);

    if (opt_index >= argc) {
        fprintf(stderr, "%s : %s\n",
                        PRG_NAME, "[err] specify a record string to process!");
        exit(1);
    }

    if (argv[1])
      strncpy (record, argv[opt_index], MAX_CHAR);

    if ( strlen(record) == 0 ) {
        fprintf(stdout, "Please enter a valid string to search pattern for!\n");
        fprintf(stdout, "Usage: %s %s %s %s\n", PRG_NAME, "[COSTS]", "[PATTERN]", "[STRING]");
        exit(1);
    }

/* Testing patterns */
//    char record[] = "ATTTACTATGTAAAGATAGAAGGAATAAGGTGAAG";
//    char regexp[]   = "GATTT";
//    char regexp[]   = "GTAAAGA";
//    char regexp[]   = "(G|T)A*GAT";

    int comp_flags  = REG_EXTENDED | REG_ICASE ;

    static regex_t preg;	      /* Compiled pattern to search for. */
    static regaparams_t match_params; /* regexp matching parameters */
    int errcode;

    regmatch_t pmatch = { 0, 0 };     /* matched pattern structure */
    regamatch_t match;                /* overall match structure */

    memset(&match, 0, sizeof(match)); /* initialize the overall match struct */
    match.pmatch = &pmatch;           /* assign default pattern structure */
    match.nmatch = 1;                 /* initialization of pmatch array */

    /* setup the default match parameters */
    tre_regaparams_default(&match_params);
    /* Set the maximum number of errors allowed for a record to match. */
    match_params.max_cost   = max_cost;
    match_params.cost_ins   = cost_ins;
    match_params.cost_del   = cost_del;
    match_params.cost_subst = cost_subs;

    fprintf(stdout, "Default regex params set by TRE:\n");
    fprintf(stdout, "\t%-12s : %4d\n", "cost_ins", match_params.cost_ins);
    fprintf(stdout, "\t%-12s : %4d\n", "cost_del", match_params.cost_del);
    fprintf(stdout, "\t%-12s : %4d\n", "cost_substr", match_params.cost_subst);
    fprintf(stdout, "\t%-12s : %4d\n", "max_cost", match_params.max_cost);
    fprintf(stdout, "\t%-12s : %4d\n", "max_ins", match_params.max_ins);
    fprintf(stdout, "\t%-12s : %4d\n", "max_del", match_params.max_del);
    fprintf(stdout, "\t%-12s : %4d\n", "max_subst", match_params.max_subst);
    fprintf(stdout, "\t%-12s : %4d\n", "max_err", match_params.max_err);
    fprintf(stdout, "\n\n");

    /* Step 1: compile the regex */
    errcode = tre_regcomp(&preg, regexp, comp_flags);
    if (errcode)
      {
        char errbuf[256];
        tre_regerror(errcode, &preg, errbuf, sizeof(errbuf));
        fprintf(stderr, "%s: %s: %s\n",
  	      PRG_NAME, "Error in search pattern", errbuf);
        exit(1);
      }

//    if (tre_regexec(&delim, "", 0, NULL, 0) == REG_OK)
//    {
//      fprintf(stderr, "%s: %s\n", PRG_NAME,
//	      "Record delimiter pattern must not match an empty string");
//      exit(1);
//    }

    /* Step 2: search for the pattern in the haystack */

    errcode = tre_regaexec(&preg, record, &match, match_params, 0);
    if (errcode == REG_OK) {
        fprintf(stdout, "Found match!\n");
        fprintf(stdout, "\t%10s : %s\n", "record", record);
        fprintf(stdout, "\t%10s : %s\n", "pattern", regexp);
        fprintf(stdout, "\t%10s : %4d\n", "cost", match.cost);
        fprintf(stdout, "\t%10s : %4d\n", "num_ins", match.num_ins);
        fprintf(stdout, "\t%10s : %4d\n", "num_del", match.num_del);
        fprintf(stdout, "\t%10s : %4d\n", "num_subst", match.num_subst);
        fprintf(stdout, "\t%10s : %2d - %2d\n", "char pos", 
                                                pmatch.rm_so, pmatch.rm_eo);
    }
    else {
        fprintf(stdout, "Found no matches!\n");
        fprintf(stdout, "%6s : %s\n", "regexp", regexp);
        fprintf(stdout, "%6s : %s\n", "record", record);

    }

    return 0;
}

int
process_options(int argc, char *argv[]) {
    int c;
    int index;
    char *opt_p_value = NULL;

    while( (c = getopt(argc, argv, "hp:D:I:S:m:")) != -1 ) {
        switch(c) {
            case 'h':
                help_message();
                exit(0);
                break;
            case 'p':
                opt_p_value = optarg;
                break;
            case 'm':
                max_cost = atoi(optarg);
                break;
            case 'I':
                cost_ins = atoi(optarg);
                break;
            case 'D':
                cost_del = atoi(optarg);
                break;
            case 'S':
                cost_subs = atoi(optarg);
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
        exit(1);
    }
    else {
        strncpy(regexp, opt_p_value, MAX_CHAR);
    }

    return optind;
}

void
help_message() {
    fprintf(stdout, "Usage: %s %s %s %s\n",
                    PRG_NAME, "[options]", "-p <pattern>", "<STRING>");
    fprintf(stdout, "\t%-20s%-20s\n", "-h", "This help message");
    fprintf(stdout, "\t%-20s%-20s\n", "-p", "Pattern of interest to grep [REQUIRED]");
    fprintf(stdout, "\t%-20s%-20s\n", "-m", "Total number of mismatches to at most allow for");
    fprintf(stdout, "\t%-20s%-20s\n", "", "in search pattern [Default: 0]");
    fprintf(stdout, "\t%-20s%-20s\n", "-I", "Cost of base insertions in obtaining");
    fprintf(stdout, "\t%-20s%-20s\n", "", "approximate match [Default: 1]");
    fprintf(stdout, "\t%-20s%-20s\n", "-D", "Cost of base deletions in obtaining");
    fprintf(stdout, "\t%-20s%-20s\n", "", "approximate match [Default: 1]");
    fprintf(stdout, "\t%-20s%-20s\n", "-S", "Cost of base substitutions in obtaining");
}
