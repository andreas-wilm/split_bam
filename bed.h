/* Trimmed down version of Devon Ryan's original https://github.com/dpryan79/MethylDackel/blob/master/MethylDackel.h
 * commit c04fb45
 */
#include <inttypes.h>
#include <stdio.h>
#include <zlib.h>
#include "htslib/sam.h"
//#include "htslib/faidx.h"


/*! @typedef
 @abstract Structure to hold one region defined in a BED file
 @field	tid	The chromosome ID, defined by bam_hdr_t
 @field start	0-based start position
 @field end	1-based end position
 @field	strand	0: Ignore strand (".")
                1: Top strand ("+")
                2: Bottom strand ("-")
 @discussion The start and end coordinates will be truncated if they go beyond
 the bounds of a chromosome, as indicated by bam_hdr_t.
*/
typedef struct {
    int32_t tid, start, end; //0-based, [start, end)
    int8_t strand; //0: ., 1: +, 2: -
} bedRegion;

/*! @typedef
 @abstract Structure to hold one region defined in a BED file
 @field region	Pointer to the regions
 @field	n	Number of regions
 @field	m	maxmum number of regions that can currently be stored.
 @discussion You must free this with destroyBED()
*/
typedef struct {
    bedRegion *region;
    int32_t n, m; //Current (n) and maximal possible (m) number of regions
} bedRegions;


//bed.c
int posOverlapsBED(int32_t tid, int32_t pos, bedRegions *regions, int idxBED);
int spanOverlapsBED(int32_t tid, int32_t start, int32_t end, bedRegions *regions, int *idx);
void sortBED(bedRegions *regions);
void destroyBED(bedRegions *regions);
bedRegions *parseBED(char *fn, bam_hdr_t *hdr);
