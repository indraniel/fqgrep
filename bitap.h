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
#define PRG_NAME "fqgrep"
#define PATTERN_LENGTH_LIMIT 64

/* P R O T O T Y P E S *******************************************************/
const char* bitap_fuzzy_bitwise_search(const char *text, 
                                        const char *pattern, 
                                        int k);

#ifdef __cplusplus
}
#endif

#endif /* _BM_H */
