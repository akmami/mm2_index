TARGET := segram
SRCS := $(wildcard *.c)
OBJS := $(SRCS:.c=.o)
CURRENT_DIR := $(shell pwd)

CC := cc
CFLAGS := -g -Wall -O2 -Wc++-compat -DHAVE_KALLOC
LDFLAGS := -lm -pthread -lz

# ========================================
#  Build Rules
# ========================================

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(CXXLIBS) $(LDFLAGS)
	rm -f $(OBJS)

%.o: %.c
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

# ========================================
#  Utility Targets
# ========================================

clean: 
	@echo "Cleaning"
	rm -f $(OBJS)
	rm -f $(TARGET)

install: clean $(TARGET)

# vg construct -r human_v38.chr22.fasta -v hprc-v1.0-pggb.grch38.1-22+X.chr22.vcf > hg38.hprc.chr22.vg
# vg ids -s hg38.hprc.chr22.vg > hg38.hprc.chr22.sorted.vg
# vg view -g hg38.hprc.chr22.sorted.vg | vg mod -O - > hg38.hprc.chr22.sorted.gfa