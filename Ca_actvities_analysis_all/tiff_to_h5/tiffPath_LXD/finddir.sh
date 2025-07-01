#!/bin/bash

find . -iregex  ".*r\s*a\s*w\s*t\s*i\s*f\s*f\s*$" -type d> targetlist.csv # for rawtiff

# find . -iregex  ".*r\s*a\s*w\s*t\s*i\s*f\s*f\s*/z[0-9]+$" -type d > targetlist.csv# for rawtiff/z01....
