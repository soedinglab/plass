#include "read.h"



//note: N->N, S->S, W->W, U->A, T->A
static const char* complement_tab =
        "................................................................"
                ".TVGH..CD..M.KN...YSAABW.R.......tvgh..cd..m.kn...ysaabw.r......"
                "................................................................"
                "................................................................";

/* Complements a canonical ASCII base (A, C, G, T, N). */
static inline char
complement(char c)
{

    return complement_tab[(unsigned char)c];
}


static inline char
identity_mapping(char c)
{
    return c;
}

/* Reverse a sequence of @len chars, applying the mapping function @map_char to
 * each, including the middle char if @len is odd.  */
static inline void
reverse_with_mapping(char *p, size_t len, char (*map_char)(char))
{
    char tmp;
    char *pp = p + len;
    while (pp > p) {
        --pp;
        tmp = *p;
        *p = (*map_char)(*pp);
        *pp = (*map_char)(tmp);
        ++p;
    }
}

/* Reverse-complement a read in place.  */
extern void
reverse_complement(struct read *r)
{
    reverse_with_mapping(r->seq, r->seq_len, complement);
    reverse_with_mapping(r->qual, r->seq_len, identity_mapping);
}
