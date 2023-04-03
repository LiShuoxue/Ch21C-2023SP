mkdir anneal
cd anneal

filepath="https://raw.githubusercontent.com/LiShuoxue/Ch21C-2023SP/main/Gromacs/Argon-gas/"

curl -O "${filepath}argon.top"
curl -O "${filepath}argon_start.pdb"
curl -O "${filepath}md.mdp" | sed "s/MYTEMP/$1/g" > mdwork.mdp
curl -O "${filepath}anneal.mdp" | sed "s/STARTTEMP/$1/g" | sed "s/ENDTEMP/$2/g" > anneal.mdp

gmx grompp -f mdwork.mdp -c argon_start.pdb -p argon.top
gmx mdrun -s topol.tpr -v -c argon_1ns.gro -nice 0
mv traj_comp.xtc gas.xtc
gmx grompp -f anneal.mdp -c argon_start.pdb -p argon.top -maxwarn 1

# Post-processing
echo 4  | gmx energy -o Potential.xvg
echo 5  | gmx energy -o Kinetic.xvg
echo 6  | gmx energy -o Total_E.xvg
echo 8  | gmx energy -o Temperature.xvg
echo 10 | gmx energy -o Pressure.xvg

cd ..
