/*  
Author: Andreas Wilm <wilma@gis.a-star.edu.sg>
License: MIT

Based on htslib test/test_view.c.
(cram removed, because of compilation problems with htslib-1.3)

Original copyright:
    Copyright (C) 2012 Broad Institute.
    Copyright (C) 2013-2014 Genome Research Ltd.
    Author: Heng Li <lh3@sanger.ac.uk>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.  */

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
                  
#include "htslib/sam.h"
#include "bed.h"

#define MYNAME "split_bam_by_chr"

                  
void usage()
{
  fprintf(stderr, "Usage: %s [-bSI] [-L compr_lvl] [-Z hdr_nuls] -l bedfile -o outprefix <in.[bam|sam]>\n", MYNAME);
  fprintf(stderr, "Options: \n");
  fprintf(stderr, "  -b : output is BAM\n");
  fprintf(stderr, "  -S : input is SAM\n");
  fprintf(stderr, "  -I : ignore SAM error\n");
  fprintf(stderr, "  -L : BAM compression level (0-9) (forces BAM output)\n");
  fprintf(stderr, "  -Z : extra header text\n");
  fprintf(stderr, "  -l : BED file listing regions (zero offset, half-open)\n");
  fprintf(stderr, "  -o : output prefix (final: prefix.region.bam)\n");
}


int main(int argc, char *argv[])
{
    int i=0;
    samFile *in;
    char *bedfile = NULL;
    bedRegions *bedregions = NULL;
    
    int flag = 0, c, clevel = -1, ignore_sam_err = 0;
    char moder[8];
    bam_hdr_t *hdr;
    bam1_t *b;
    htsFile **out;
    char *outprefix = 0;
    char modew[800];
    int r = 0, exit_code = 0;
    int extra_hdr_nuls = 0;
    int unaln_idx = -1;
    while ((c = getopt(argc, argv, "bSIl:L:Z:o:")) >= 0) {
      switch (c) {
        case 'b': flag |= 2; break;/* bam out */
        case 'S': flag |= 1; break;/* sam in */
        case 'I': ignore_sam_err = 1; break;
        case 'l': bedfile = strdup(optarg); break;
        case 'L': clevel = atoi(optarg); flag |= 2; break;/* compression level */
        case 'Z': extra_hdr_nuls = atoi(optarg); break;
        case 'o': outprefix = strdup(optarg); break;
        }
    }
    if (argc == optind) {
      usage();
      return 1;
    }

    if (! outprefix) {
      fprintf(stderr, "Need output prefix\n");
        return EXIT_FAILURE;
    }
    if (! bedfile) {
      fprintf(stderr, "Need bed file\n");
        return EXIT_FAILURE;
    }

    /* setting up input
     */
    strcpy(moder, "r");
    if (flag&4) strcat(moder, "c");
    else if ((flag&1) == 0) strcat(moder, "b");

    in = sam_open(argv[optind], moder);
    if (in == NULL) {
        fprintf(stderr, "Error opening \"%s\"\n", argv[optind]);
        return EXIT_FAILURE;
    }
    hdr = sam_hdr_read(in);
    if (hdr == NULL) {
        fprintf(stderr, "Couldn't read header for \"%s\"\n", argv[optind]);
        return EXIT_FAILURE;
    }
    hdr->ignore_sam_err = ignore_sam_err;
    if (extra_hdr_nuls) {
        char *new_text = realloc(hdr->text, hdr->l_text + extra_hdr_nuls);
        if (new_text == NULL) {
            fprintf(stderr, "Error reallocing header text\n");
            return EXIT_FAILURE;
        }
        hdr->text = new_text;
        memset(&hdr->text[hdr->l_text], 0, extra_hdr_nuls);
        hdr->l_text += extra_hdr_nuls;
    }
    b = bam_init1();

    /* read bed file
     */
    assert(bedfile);
    bedregions = parseBED(bedfile, hdr);
    if (! bedregions) {
      fprintf(stderr, "Error parsing of %s failed\n", bedfile);
      return EXIT_FAILURE;
    }
    
    /* setting up output
     */
    strcpy(modew, "w");
    if (clevel >= 0 && clevel <= 9) sprintf(modew + 1, "%d", clevel);
    if (flag&8) strcat(modew, "c");
    else if (flag&2) strcat(modew, "b");
    unaln_idx = bedregions->n;
    out = calloc(bedregions->n+1, sizeof(htsFile*));/* last one is for unaligned reads */
    for (i=0; i<bedregions->n+1; i++) {/* +1 incl unaliged */
      char *buffer = malloc(1024 * sizeof(char));
      assert(buffer);
      if (i==unaln_idx) {
        snprintf(buffer, 1024, "%s.unaligned.bam", outprefix);
      } else {
        const char *sq = hdr->target_name[bedregions->region[i].tid];
        snprintf(buffer, 1024, "%s.%s:%" PRId32 "-%" PRId32 ".bam", outprefix, sq,
                 bedregions->region[i].start, bedregions->region[i].end);
      }
#ifdef DEBUG
      fprintf(stderr, "DEBUG: writing to %s\n", buffer);
#endif
      out[i] = hts_open(buffer, modew);
      if (out[i] == NULL) {
        fprintf(stderr, "Error opening %s\n", buffer);
        return EXIT_FAILURE;
      }
      free(buffer);
    }

    /* writing 
     */
    for (i=0; i<bedregions->n+1; i++) {
      if (sam_hdr_write(out[i], hdr)) {
        fprintf(stderr, "Error sam header writing failed.\n");
        exit(1);
      }
    }
    while ((r = sam_read1(in, hdr, b)) >= 0) {
      int write_error = 0;

      if (b->core.flag & BAM_FUNMAP) {
            if (sam_write1(out[unaln_idx], hdr, b) < 0) {
              fprintf(stderr, "Error writing output.\n");
              write_error = exit_code = 1;
            }        
      } else {
        /* naive looping over regions */
        for (i=0; i<bedregions->n; i++) {
          int ovlp = spanOverlapsBED(b->core.tid, b->core.pos, bam_endpos(b), bedregions, &i);
          if (ovlp==1) {
            if (sam_write1(out[i], hdr, b) < 0) {
              fprintf(stderr, "Error writing output.\n");
              write_error = exit_code = 1;
            }
            break;/* since regions should be exclusive we can exit here */
          }
        }
      }
      if (write_error) break;
    }
    if (r < -1) {
        fprintf(stderr, "Error parsing input.\n");
        exit_code = 1;
    }

    for (i=0; i<bedregions->n+1; i++) {
      r = sam_close(out[i]);
      if (r < 0) {
        fprintf(stderr, "Error closing output.\n");
        exit_code = 1;
      }
    }

    bam_destroy1(b);
    bam_hdr_destroy(hdr);
    destroyBED(bedregions);
    r = sam_close(in);
    if (r < 0) {
        fprintf(stderr, "Error closing input.\n");
        exit_code = 1;
    }

    free(out);
    free(outprefix);
    free(bedfile);

    return exit_code;
}
