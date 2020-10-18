#!/bin/bash

input_cov_file=$1;
i=0;

cat $input_cov_file | while read line
do
	let i+=1;
	a=$(echo "$line" | sed 's/[:|]\+/\t/g' | awk '{ if( $17 != 0) { print $10/$17 } }')
	b=$(echo "$line" | sed 's/[:|]\+/\t/g' | awk '{ if( $17 != 0) { print $12/$17 } }')
	c=$(echo "$line" | sed 's/[:|]\+/\t/g' | awk '{ if( $17 != 0) { print $14/$17 } }')
	d=$(echo "$line" | sed 's/[:|]\+/\t/g' | awk '{ if( $17 != 0) { print $16/$17 } }')
	maximum=$(echo -e "$a\n$b\n$c\n$d" | sort -rn | head -1)
	echo $i $line $maximum
done

