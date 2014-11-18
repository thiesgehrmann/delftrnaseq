#!/bin/sh

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

IS_mean=0;
IS_stdev=0;
for s in "$@"; do
  dist=`head -n 100000 $s | grep -e '^[^@]' | awk 'and(2, $2) && and(128, $2)' | cut -f9 | head -n20000 | "$SCRIPT_DIR"/stats.sh`
  c_IS_mean=`echo $dist | cut -d\  -f1`;
  c_IS_stdev=`echo $dist | cut -d\  -f2`;
  IS_mean=`echo $IS_mean + $c_IS_mean | bc`;
  IS_stdev=`echo $IS_stdev + $c_IS_stdev | bc`;
done
IS_mean=`echo $IS_mean / $# | bc`;
IS_stdev=`echo $IS_stdev / $# | bc`;

echo -en "$IS_mean $IS_stdev"
