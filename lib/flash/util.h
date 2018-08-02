#ifndef _FLASH_UTIL_H_
#define _FLASH_UTIL_H_

#include <stddef.h>
#include <pthread.h>
#include <stdio.h>

#define ARRAY_LEN(A) (sizeof(A) / sizeof((A)[0]))

#ifdef __GNUC__
#	if __GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 4)
#		define __cold __attribute__((cold))
#	else
#		define __cold
#	endif
#	define __noreturn __attribute__((noreturn))
#	define __format(type, format_str, args_start) \
			__attribute__((format(type, format_str, args_start)))
#	define max(a,b) ({ __typeof__(a) _a = (a); __typeof__(b) _b = (b); _a > _b ? _a : _b; })
#	define min(a,b) ({ __typeof__(a) _a = (a); __typeof__(b) _b = (b); _a < _b ? _a : _b; })
#	define inline inline __attribute__((always_inline))
#else
#	define __noreturn
#	define __cold
#	define __format(type, format_str, args_start)
#	define max(a,b) (((a) > (b)) ? (a) : (b))
#	define min(a,b) (((a) < (b)) ? (a) : (b))
#endif


extern void
fatal_error(const char *msg, ...) __noreturn __cold __format(printf, 1, 2);

extern void
fatal_error_with_errno(const char *msg, ...) __noreturn __cold __format(printf, 1, 2);

extern unsigned long warning_count;

extern void
warning(const char *msg, ...) __cold __format(printf, 1, 2);

extern FILE *infofile;

extern void
info(const char *msg, ...) __format(printf, 1, 2);

extern void *
xmalloc(size_t size);

#ifdef NDEBUG
#  define xfree(p, size) free(p)
#else
extern void
xfree(void *p, size_t size);
#endif

extern void *
xzalloc(size_t size);

extern char *
xstrdup(const char *str);

extern void *
xrealloc(void *ptr, size_t size);

extern unsigned
get_default_num_threads(void);

extern void
mkdir_p(const char *dir);

extern pthread_t
create_thread(void *(*proc)(void *), void *params);

extern void
join_thread(pthread_t t);


#endif /* _FLASH_UTIL_H_ */
