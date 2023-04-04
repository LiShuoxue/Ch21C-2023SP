import os, subprocess
import pandas as pd
import numpy as np

units = {
    "Potential": "kJ/mol",
    "Kinetic": "kJ/mol",
    "Total_E": "kJ/mol",
    "Temperature": "K",
    "Pressure": "Bar"
}

def generate_Argon_pdbfile(a, N, filename=None):

    if filename is None: filename = "argon-start-N{}-a{:.3f}.pdb".format(N, a)

    with open(filename+".pdb", "w") as f:
        f.write("HEADER    Argon\n")
        f.write("REMARK    Box with a = {:.3f}, N = {}\n".format(a, N))
        f.write("CRYST1   {:.3f}   {:.3f}   {:.3f}  90.00  90.00  90.00 P 1  1\n".format(a,a,a))
        f.write("MODEL        1\n")
        for i in range(N):
            f.write("ATOM\t{}\tAr\tAr\t{}\t{:.3f}\t{:.3f}\t{:.3f}\t1.00\t0.00\n".format(i+1, i+1, np.random.random() * a, np.random.random() * a, np.random.random() * a))
        f.write("TER\n")
        f.write("ENDMDL\n")
        f.close()

def run_mds(tmps):

    with open("runmd.sh", 'w') as f:
        f.write('''# Shuoxue Li <sli7@caltech.edu>
# Download sample file from github
filepath="https://raw.githubusercontent.com/LiShuoxue/Ch21C-2023SP/main/Gromacs/Argon-gas/"
for filename in "argon_start.pdb" "argon.top" 
do
  curl -O ${filepath}${filename}
done

# change temperature to the arg $1
curl -O ${filepath}md.mdp
cat md.mdp | sed "s/MYTEMP/$1/g" > mdwork.mdp

# Run Molecular dynamics
gmx grompp -f mdwork.mdp -c argon_start.pdb -p argon.top
gmx mdrun -s topol.tpr -v -c argon_1ns.gro -nice 0

# Post-processing
echo 4  | gmx energy -o Potential.xvg
echo 5  | gmx energy -o Kinetic.xvg
echo 6  | gmx energy -o Total_E.xvg
echo 8  | gmx energy -o Temperature.xvg
echo 10 | gmx energy -o Pressure.xvg
''')
    
    df = pd.DataFrame()
    for tag in ['Potential', 'Kinetic', 'Total_E', 'Temperature', 'Pressure']:
        df[tag] = np.zeros(len(tmps))
    for i, tmp in enumerate(tmps):

        folder = "T-{:.2f}K".format(tmp)

        try: os.makedirs(folder)
        except FileExistsError: pass
        os.chdir(folder)

        subprocess.call(["sh", "../runmd.sh", "{:.2f}".format(tmp)])
        for tag in ['Potential', 'Kinetic', 'Total_E', 'Temperature', 'Pressure']:
            step, quant = np.loadtxt("{}.xvg".format(tag), comments=['#', '@']).T
            quant_mean = np.mean(quant)
            quant_std  = np.std(quant)

            print("Point {} : {} = ({} ± {}) {}".format(
                i+1, tag, quant_mean, quant_std, units[tag]
            ))

            df[tag][i] = quant_mean

        os.chdir("..")
        print("-" * 50)

    df.to_csv("results.csv")
