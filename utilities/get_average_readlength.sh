#!/bin/sh

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )";


alignlen=`cat "$@" | grep 'Average mapped length' | cut -d\| -f2 | sed -e 's/[^0-9.]//' | "$SCRIPT_DIR"/stats.sh | cut -d\  -f1`;
read_length=`echo $alignlen / 2 | bc`;

echo -en "$read_length"

