#!/bin/bash

. /home/mario/software/heasoft-6.26.1/x86_64-pc-linux-gnu-libc2.27/headas-init.sh
echo "HEASoft has been initiated."

if [ $# -lt 2 ]; then
    echo "$0: Not enough arguments. Inputs are 'original gti' and 'new ascii gti data'."
    exit 2
elif [ $# -gt 2 ]; then
    echo "$0: Too many arguments. Arguments are 'original gti' and 'new ascii gti data'."
    exit 2
fi

old_gti=$1
new_gti=$2

flcol $old_gti outfile=lcusergti_cols.txt clobber=yes
wait
fcreate cdfile=lcusergti_cols.txt datafile=$new_gti extname="USER_GTI" outfile=lcusergti.fits clobber=yes

echo "'lcusergti.fits' created."

