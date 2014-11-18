#!/bin/sh

awk '{ count += 1; sum += $1; sumsq+=$1*$1 } END { print sum/count, sqrt(sumsq/count - (sum/count)**2) }'

