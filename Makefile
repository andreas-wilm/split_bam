 # user args
CC?=gcc
HTSLIBDIR?=/mnt/projects/rpd/apps.testing/htslib-1.3/;
#HTSLIBDIR?=/Users/wilma/local/dotkit-inst/htslib-1.3/

CFLAGS = -Wall -O3 -DNDEBUG $(IFLAGS) -I$(HTSLIBDIR)/include/  -I$(HTSLIBDIR)/
LFLAGS = -lz -pthread -lcrypto
LFLAGS_STATIC = $(LFLAGS)
LFLAGS_DYN = $(LFLAGS) -lhts -L$(HTSLIBDIR)/
LIBHTSA = $(HTSLIBDIR)/libhts.a
static: split_bam.c $(LIBHTSA)
	$(CC) $(CFLAGS) split_bam.c bed.c -o split_bam $($IFLAGS) $(LIBHTSA) $(LFLAGS_STATIC) -lcurl 

dynamic: split_bam.c
	$(CC) $(CFLAGS) split_bam.c bed.c -o split_bam $($IFLAGS) $(LFLAGS_DYN)

all: static
