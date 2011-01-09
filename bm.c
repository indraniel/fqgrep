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

/* I N C L U D E S ***********************************************************/
# include "bm.h"

/* F U N C T I O N S *********************************************************/
void 
compute_prefix(const char* str, size_t size, int result[size]) {
    size_t q;
    int k;
    result[0] = 0;

    k = 0;
    for (q = 1; q < size; q++) {
        while (k > 0 && str[k] != str[q])
            k = result[k-1];

        if (str[k] == str[q])
            k++;
        result[q] = k;
    }
}
 
void 
prepare_badcharacter_heuristic(const char *str, 
                               size_t size,
                               int result[ALPHABET_SIZE]) {
 
    size_t i;

    for (i = 0; i < ALPHABET_SIZE; i++)
        result[i] = -1;

    for (i = 0; i < size; i++)
        result[(size_t) str[i]] = i;
}
 
void prepare_goodsuffix_heuristic(const char *normal, 
                                  size_t size,
                                  int result[size + 1]) {
    char *left = (char *) normal;
    char *right = left + size;
    char reversed[size+1];
    char *tmp = reversed + size;
    size_t i;

    /* reverse string */
    *tmp = 0;
    while (left < right)
        *(--tmp) = *(left++);

    int prefix_normal[size];
    int prefix_reversed[size];

    compute_prefix(normal, size, prefix_normal);
    compute_prefix(reversed, size, prefix_reversed);

    for (i = 0; i <= size; i++) {
        result[i] = size - prefix_normal[size-1];
    }

    for (i = 0; i < size; i++) {
        const int j = size - prefix_reversed[i];
        const int k = i - prefix_reversed[i]+1;

        if (result[j] > k)
            result[j] = k;
    }
}

const char*
boyermoore_search(const char *haystack, const char *needle) {
    /*
    * Calc string sizes
    */
    size_t needle_len, haystack_len;
    needle_len = strlen(needle);
    haystack_len = strlen(haystack);

    /*
    * Simple checks
    */
    if(haystack_len == 0)
        return NULL;
    if(needle_len == 0)
        return haystack;

    /*
    * Initialize heuristics
    */
    int badcharacter[ALPHABET_SIZE];
    int goodsuffix[needle_len+1];

    prepare_badcharacter_heuristic(needle, needle_len, badcharacter);
    prepare_goodsuffix_heuristic(needle, needle_len, goodsuffix);

    /*
    * Boyer-Moore search
    */
    size_t s = 0;
    while(s <= (haystack_len - needle_len))
    {
        size_t j = needle_len;
        while(j > 0 && needle[j-1] == haystack[s+j-1])
            j--;

        if(j > 0)
        {
            int k = badcharacter[(size_t) haystack[s+j-1]];
            int m;
            if(k < (int)j && (m = j-k-1) > goodsuffix[j])
                s+= m;
            else
                s+= goodsuffix[j];
        }
        else
        {
            return haystack + s;
        }
    }

    /* not found */
    return NULL;
}
