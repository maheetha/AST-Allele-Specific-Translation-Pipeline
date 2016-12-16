#!/usr/bin/bash

<<COMMENTS

The purpose is to jsut get the first 12 characters of the sample file and add it. 

COMMENTS
cut -c1-12 $1 > $1.edited

cut -f3- $1 > $1.edited2

paste -d'\t' $1.edited $1.edited2  > $1.finaledited


