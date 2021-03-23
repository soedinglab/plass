#include "read.h"

//note: N->N, S->S, W->W, U->A, T->A
static const char* complement_tab =
        "................................................................"
        ".TVGH..CD..M.KN...YSAABW.R.......tvgh..cd..m.kn...ysaabw.r......"
        "................................................................"
        "................................................................";

/* Complements a canonical ASCII base (A, C, G, T, N). */
struct complement {
    char operator()(char c) {
        return complement_tab[(unsigned char)c];
    }
};


struct identity_mapping {
    char operator()(char c) {
        return c;
    }
};

/* Reverse a sequence of @len chars, applying the mapping function @map_char to
 * each, including the middle char if @len is odd.  */
template<typename T>
static inline void
reverse_with_mapping(char *p, size_t len)
{
    T t;
    char tmp;
    char *pp = p + len;
    while (pp > p) {
        --pp;
        tmp = *p;
        *p = t(*pp);
        *pp = t(tmp);
        ++p;
    }
}

/* Reverse-complement a read in place.  */
extern void
reverse_complement(struct read *r)
{
    reverse_with_mapping<complement>(r->seq, r->seq_len);
    reverse_with_mapping<identity_mapping>(r->qual, r->seq_len);
}
