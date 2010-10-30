/* N O T E S *****************************************************************/
/* compliation line:
   gcc -g -o tre-test -ltre tre-test.c
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

/* M A I N *******************************************************************/
int main(int argc, char *argv[]) {
    char record[MAX_CHAR] = "";
    char regexp[MAX_CHAR] = "";

    if (argv[1])
      strncpy (regexp, argv[1], MAX_CHAR);

    if ( strlen(regexp) == 0 ) {
        fprintf(stdout, "Please enter a regexp pattern!\n");
        fprintf(stdout, "Usage: %s %s %s\n", PRG_NAME, "[PATTERN]", "[STRING]");
        exit(1);
    }

    if (argv[2])
      strncpy (record, argv[2], MAX_CHAR);

    if ( strlen(record) == 0 ) {
        fprintf(stdout, "Please enter a string to search pattern for!\n");
        fprintf(stdout, "Usage: %s %s %s\n", PRG_NAME, "[PATTERN]", "[STRING]");
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
    match_params.max_cost = 1;

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
