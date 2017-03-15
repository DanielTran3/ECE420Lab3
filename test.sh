#! /bin/bash
# Run your client for 100 times
# 
# 
# Usage:
#	./test.sh
# Notes:
#	This shell script should be in the same directory as your
#	client implementation and your client should have portnumber 
#	and the number of strings in theArray as command line paramters;
#


# Parameters
Duplicates=50


clear

echo "Start..."
ATTEMPT=0
while [[ $ATTEMPT -ne $Duplicates ]]; do
	let ATTEMPT+=1
	./gje_seq
	sleep .5
done
