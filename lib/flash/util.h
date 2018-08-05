#ifndef _FLASH_UTIL_H_
#define _FLASH_UTIL_H_

#include <stddef.h>
#include <pthread.h>
#include <stdio.h>

#define ARRAY_LEN(A) (sizeof(A) / sizeof((A)[0]))

extern void *
xmalloc(size_t size);

extern void *
xrealloc(void *ptr, size_t size);

#endif /* _FLASH_UTIL_H_ */
