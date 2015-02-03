#!/bin/sh
for Dir in $(find processor* -type d ); 
do
    cp streams.in mech.xml $Dir
done
