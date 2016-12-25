#!/bin/bash
echo "Bash version ${BASH_VERSION}..."
for i in {1..9..1}
  do
SEARCH="current_a"
OUTPUT="deneme"

cat temp_Parameters | sed -e "s/$SEARCH/$i/" > Parameters.h #changes current_a word in temp_Parameter with $i and writes whole text to the Parameters.h 
cat temp_Makefile | sed -e "s/$OUTPUT/"a="$i/" > Makefile #changes deneme word in temp_Makefile with a=$i and writes whole text to the Makefile
make #creates the executable file for relevant a value.
gnome-terminal -e "bash -c \"./a=$i; exec bash\"" #opens new terminal and executes the line ./a=$i there and continues to the next value in the loop 
 done


