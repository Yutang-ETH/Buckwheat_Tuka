#!/bin/bash

yak count -k31 -b37 -t24 -o pat.yak <(zcat 20230213.X-FE228_R1.fastq.gz) <(zcat 20230213.X-FE228_R2.fastq.gz)
yak count -k31 -b37 -t24 -o mat.yak <(zcat 20230213.X-FE213_1_R1.fastq.gz) <(zcat 20230213.X-FE213_1_R2.fastq.gz)
