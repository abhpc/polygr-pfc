#!/bin/bash -e
density=$1
width=$2
height=$3
project_dir=$PFC_ROOT
coordinates="${density%.dat}.xy"
neighbors="${density%.dat}.nh"

source /etc/bashrc
module purge
module load openjdk/13.0.1
java -jar "$project_dir/java/coordinator.jar" "$density" "$width" "$height" 0.7 0.7 7.3 2.46 "$coordinates" "$neighbors" > conversion.stdout 2> conversion.stderr
awk 'END {print NR}' "$coordinates" > atoms.count
{
	cat atoms.count
	printf '\n'
	awk '{print $1+1,1,$2,$3,0}' "$coordinates"
} > graphene.xyz

module purge
module load ovito/2.9.0
ovitos "$project_dir/source/xyz2lmp.py" >> conversion.stdout 2>> conversion.stderr
sed -i 's@0.0 0.0 zlo zhi@-10.0 10.0 zlo zhi@g' graphene.data
