#! bin/bash

#$ -j y
#$ -cwd
#$ -m e
#$ -q "short.q"

./myQuadTree -p 4096 -d 4 -i idx4 -f ../DataGen/RectUni4D100K.txt -q Queries100K.txt -r 20 -m AA -h 7 -t 0 -o 0 > OutputUni4D100K.txt

