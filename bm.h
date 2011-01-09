/* L I C E N S E *************************************************************/
/*

   Copyright (C) 2010, 2011 Wikipedia®
   http://en.wikipedia.org/wiki/Wikipedia:About

   This work is licensed under a Creative Commons Attribution-ShareAlike
   3.0 Unported License (http://creativecommons.org/licenses/by-sa/3.0/)
*/

/* N O T E S *****************************************************************/
/* Boyer-Moore search algorithm

   Taken from wikipedia:
   http://en.wikipedia.org/wiki/Boyer%E2%80%93Moore_string_search_algorithm
*/

#ifndef _BM_H_
#define _BM_H_

#ifdef __cplusplus
extern "C" {
#endif

/* I N C L U D E S ***********************************************************/
#include <stdlib.h>
#include <string.h>
#include <limits.h>

/* D E F I N E S *************************************************************/
#define ALPHABET_SIZE ( 1 << CHAR_BIT)

/* P R O T O T Y P E S *******************************************************/
void compute_prefix(const char* str, size_t size, int result[size]);
void prepare_badcharacter_heuristic(const char *str, 
                                    size_t size, 
                                    int result[ALPHABET_SIZE]);
void prepare_goodsuffix_heuristic(const char *normal, 
                                  size_t size, 
                                  int result[size + 1]);
const char* boyermoore_search(const char *haystack, const char *needle);

#ifdef __cplusplus
}
#endif

#endif /* _BM_H */
