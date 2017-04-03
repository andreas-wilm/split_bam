 # user args
CC?=gcc
HTSDIR?=/mnt/projects/rpd/apps.testing/htslib-1.3/;
#HTSDIR?=/Users/wilma/local/dotkit-inst/htslib-1.3/

CFLAGS = -Wall -O3 -DNDEBUG $(IFLAGS) -I $(HTSDIR)/include/ 
LFLAGS_STATIC = -lz -pthread
LFLAGS_DYN = -lz -pthread -lhts -L$(HTSDIR)/
LIBHTSA = $(HTSDIR)/libhts.a
static: split_bam.c $(LIBHTSA)
	$(CC) $(CFLAGS) split_bam.c bed.c -o split_bam $($IFLAGS) $(LIBHTSA) $(LFLAGS_STATIC) -lcurl 

dynamic: split_bam.c
	$(CC) $(CFLAGS) split_bam.c bed.c -o split_bam $($IFLAGS) $(LFLAGS_DYN)

all: static
