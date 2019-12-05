/*
 * util.c: miscellaneous useful functions for FLASH
 */

/*
 * Copyright (C) 2012 Tanja Magoc
 * Copyright (C) 2012, 2013, 2014 Eric Biggers
 *
 * This file is part of FLASH, a fast tool to merge overlapping paired-end
 * reads.
 *
 * FLASH is free software; you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free
 * Software Foundation; either version 3 of the License, or (at your option)
 * any later version.
 *
 * FLASH is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License
 * along with FLASH; if not, see http://www.gnu.org/licenses/.
 */

#include <errno.h>
#include <limits.h>
#include <pthread.h>
#include <stdarg.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include "util.h"

#ifdef __WIN32__
/* Get the pthread mutex declarations as a replacement for flockfile() and
     * funlockfile(). */
#  include <pthread.h>
    /* Get the GetSystemInfo() declaration as replacement for
     * sysconf(_SC_NPROCESSORS_ONLN). */
#  include <windows.h>
#endif


/* Like malloc(), but aborts if out of memory, and always returns non-NULL, even
 * if 0 bytes were requested.  */
void *
xmalloc(size_t size)
{
    void *p = malloc(size);
    if (p)
        return p;
    if (!size) {
        p = malloc(1);
        if (p)
            return p;
    }
    fprintf(stderr, "Out of memory: tried to allocate %zu bytes", size);
    return NULL;
}



/* Like realloc(), but aborts if out of memory, and always returns non-NULL,
 * even if 0 bytes were requested.  */
void *
xrealloc(void *ptr, size_t size)
{
    void *p = realloc(ptr, size);
    if (p)
        return p;
    if (!size) {
        p = malloc(1);
        if (p)
            return p;
    }
    fprintf(stderr, "Out of memory: tried to reallocate %zu bytes", size);
    return NULL;
}

