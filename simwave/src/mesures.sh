#!/bin/bash

while read bloc;do
	s=`./bin/simwave  -v --grid 100,100,100 --iter 10 --check --epsilon 1e-20 --bloc $bloc`
#	echo $bloc
done < taille
	
