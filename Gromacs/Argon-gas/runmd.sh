# Shuoxue Li <sli7@caltech.edu>

# Download sample file from github
filepath="https://raw.githubusercontent.com/LiShuoxue/Ch21C-2023SP/main/Gromacs/Argon-gas/"
for filename in "argon_start.pdb" "argon.top" "md.mdp"
do
  curl -O ${filepath}${filename}
done

# change temperature to the arg $1
cat md.mdp | sed "s/MYTEMP/$1/g" > md.mdp

# Run Molecular dynamics
gmx grompp -f md.mdp -c argon_start.pdb -p argon.top
gmx mdrun -s topol.tpr -v -c argon_1ns.gro -nice 0

# Post-processing
echo 4  | gmx energy -o Potential.xvg
echo 5  | gmx energy -o Kinetic.xvg
echo 6  | gmx energy -o Total_E.xvg
echo 8  | gmx energy -o Temperature.xvg
echo 10 | gmx energy -o Pressure.xvg
