#!/bin/bash

# This script groups files in the current directory by their file extension.
l="001 005 006 007 008 009 010 011 027 028 029 030 031 032 033 034 035 036 037 038 039 041 042 077 079 080"

for i in $l; do
    Grouper $i
done