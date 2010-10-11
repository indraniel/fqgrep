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
