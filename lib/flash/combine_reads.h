#ifndef _FLASH_COMBINE_READS_H_
#define _FLASH_COMBINE_READS_H_

#include <stdbool.h>
#include "read.h"
/* Parameters for the core algorithm of FLASH.  See the help output for more
 * information.  */
struct combine_params {
    /* --min-overlap  */
    int min_overlap;

    /* --max-overlap  */
    int max_overlap;

    /* --max-mismatch-density  */
    float max_mismatch_density;

    /* --cap-mismatch-quals  */
    bool cap_mismatch_quals;

    /* --allow-outies  */
    bool allow_outies;
};

/* Result of a call to combine_reads()  */
enum combine_status {
    /* The reads could not be combined.  */
            NOT_COMBINED = 0,

    /* The reads were combined in "innie" orientation, like the following:
     *
     * ---------->
     *     <------------
     *
     * (Note: read_2 is reverse complemented before the call to
     * combine_reads()).  */
            COMBINED_AS_INNIE,

    /* The reads were combined in "outie" orientation, like the following:
     *
     * <----------
     *     ------------>
     *
     * (Note: read_2 is reverse complemented before the call to
     * combine_reads()).  */
            COMBINED_AS_OUTIE,
};

extern enum combine_status
combine_reads(const struct read *read_1, const struct read *read_2,
              struct read *combined_read,
              const struct combine_params *params);

#endif /* _FLASH_COMBINE_READS_H_  */
