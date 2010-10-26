/* This software is in the public domain and has no warrantee */
/* 14 Apr 2005 : Fix the bug reported by Mike, by simplifying the main loop */
#ifndef LIBBITAP_H
#define LIBBITAP_H
#define LIBBITAP_MAJOR_VERSION 1
#define LIBBITAP_MINOR_VERSION 2

/* This structure contains data that is private to libbitap and should not */
/* be modified directly */
typedef struct {
  int l, n;            /* l is the number of states. n is the number of */
                       /* 'or'ing operations stored. */
  int *s, *x;          /* s shows which letter is allowed in which place. */
                           /* x determines in which areas exact matches are */
  struct bitapJumpType *j; /* required. */
  int approxMode;          /* Only used during parsing. To determine x */
} bitapType;

int NewBitap (bitapType *b, const char *regex);
const char *FindWithBitap (bitapType *b, const char *in, int inl, int e,
			   int *ereturn, const char **bReturn);
void DeleteBitap (bitapType *b);

#endif
