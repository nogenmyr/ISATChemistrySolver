#!/bin/sh
for Dir in $(find processor* -type d ); 
do
    ln -s ../mech.xml $Dir/mech.xml
#    ln -s ../streams.in  $Dir/streams.in
done
