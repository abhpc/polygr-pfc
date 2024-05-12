#!/bin/bash
#SBATCH --job-name=polygrpfc

# Patition to submit
#SBATCH --partition=MX

# The number of computing nodes
#SBATCH --nodes=4

# The number of threads per node
#SBATCH --ntasks-per-node=48

# Out and error files
#SBATCH --error=%j.err
#SBATCH --output=%j.out

# Load polygrpfc module
source /etc/bashrc
module load polygrpfc/1.1

# Unset environment variable DISPLAY for running headless OVITO
export DISPLAY=""

# Define the width of graphene sheet, unit: ang
W=1000

# Define the height of graphene sheet, unit: ang
H=1000

# Number of grains
N=100

# Initial radius of seed, unit: ang
R=2

######################################################################################
# The following lines do not require modification unless their meaning is understood #
######################################################################################

# Calculate the parameters for phase field crystal program
WW=`awk 'BEGIN{print int(1024/243*'$W')}'`
HH=`awk 'BEGIN{print int(1024/243*'$H')}'`
NN=$N
RR=$R

# copy the input templates, step1.in and step2.in, into current directory.
rm -rf *.in
cp $PFC_INP/step*.in ./

# Generate a random number, and write necessary parameters into input file step1.in
RX=`head /dev/urandom|cksum| awk '{print $1}'`
sed -i "s@xabc@$RX@g" step1.in
sed -i "s@WW@$WW@g" step1.in
sed -i "s@HH@$HH@g" step1.in
sed -i "s@NN@$NN@g" step1.in
sed -i "s@RR@$RR@g" step1.in

# Set the width and height in input file step2.in
sed -i "s@WW@$WW@g" step2.in
sed -i "s@HH@$HH@g" step2.in

# Start the initial run: step1
echo "Now run Step 1:"
mpirun polygrpfc step1
echo "Step 1 done."
echo ""


# Start relaxing run: step2
echo "Now run Step 2:"
mpirun polygrpfc step2
echo "Step 2 done."
echo ""


# Convert the file "step2-t-10000.dat" from step 2, the last step (10000 steps), into an XY file, noting that the atom IDs start from 0 here.
echo "Convert the images of Step 1 and 2:"
java -jar $PFC_JAR/coordinator.jar step2-t-10000.dat $WW $HH 0.7 0.7 7.3 2.46 \
          step2-t-10000.xy step2-t-10000.nh
java -jar $PFC_JAR/plotter.jar step1-t-# step1-t-# $WW $HH 0 1000 10000
java -jar $PFC_JAR/plotter.jar step2-t-# step2-t-# $WW $HH 0 1000 10000
echo "Images done."
echo ""

# Convert the XY file to an XYZ file, where the first column of IDs is changed from starting at 0 to starting at 1 by adding 1 to the first column. 
# The second column represents the atom type ID (uniformly set to 1).
# The third and fourth columns represent the x and y coordinates, respectively.
# The fifth column represents the z coordinate, which is uniformly defined as 0.
awk '{print $1+1,1,$2,$3,0}' step2-t-10000.xy > graphene.xyz
atomN=`tail -1 graphene.xyz |awk '{print $1}'`
sed -i "1i\ " graphene.xyz
sed -i "1i$atomN" graphene.xyz


# Call the Python interpreter of OVITO, read an XYZ file, and export it as a LAMMPS data file named "graphene.data".
rm -rf xyz2lmp.py
cat << EOF > xyz2lmp.py
from ovito.io import import_file
from ovito.io import export_file
node = import_file("graphene.xyz", columns = ["Particle Identifier", "Particle Type", "Position.X", "Position.Y", "Position.Z"])
export_file(node, "graphene.data", "lammps_data")
EOF
ovitos xyz2lmp.py
rm -rf xyz2lmp.py

# Note that in the current graphene.data file, the z-direction is defined as (0.0, 0.0), but it needs to be changed to (-10, 10).
sed -i "s@0.0 0.0 zlo zhi@-10.0 10.0 zlo zhi@g" graphene.data

mv graphene.data $SLURM_JOB_NAME.data