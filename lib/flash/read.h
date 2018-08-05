#ifndef _FLASH_READ_H_
#define _FLASH_READ_H_

#include <stddef.h>

/* In-memory representation of a DNA sequence read.  */
struct read {
	/* The "tag" that identifies the read.  */
	char *tag;

	/* The sequence of the read in ASCII characters (A, C, G, T, N).  */
	char *seq;

	/* The quality scores of the read, scaled to start at 0.  */
	char *qual;

	/* Length of the tag string.  */
	int tag_len;

	/* Length of the sequence string (number of bases in the read).  */
	int seq_len;

	/* Length of the quality string (will be equal to seq_len).  */
	int qual_len;

	/* Allocated sizes of the seq, tag, and qual buffers, respectively.  */
	size_t seq_bufsz;
	size_t tag_bufsz;
	size_t qual_bufsz;
};

struct input_stream;
#include <inttypes.h>

extern void
reverse_complement(struct read *r);


#endif /* _FLASH_READ_H_ */
