import os, subprocess, argparse
import pandas as pd
import numpy as np

units = {
    "Potential": "kJ/mol",
    "Kinetic": "kJ/mol",
    "Total_E": "kJ/mol",
    "Temperature": "K",
    "Pressure": "Bar",
    "Volume": "Angstrom^3"
}

def generate_Argon_pdbfile(V, N, filename=None):

    a = V**(1/3)

    if filename is None: filename = "argon-start-N{}-V{}.pdb".format(N, V)

    with open(filename+".pdb", "w") as f:
        f.write("HEADER    Argon\n")
        f.write("REMARK    THIS IS A SIMULATION BOX\n")
        f.write("CRYST1 {:.3f} {:.3f} {:.3f}  90.00  90.00  90.00 P 1  1\n".format(a,a,a))
        f.write("MODEL        1\n")
        for i in range(N):
            f.write("ATOM{:>7}  Ar   Ar{:>6}{:>12.3f}{:>8.3f}{:>8.3f}  1.00  0.00\n".format(i+1, i+1, np.random.random()*a, np.random.random()*a, np.random.random()*a))
        f.write("TER\n")
        f.write("ENDMDL\n")
        f.close()

def run_mds(tmps, Vs, Ns):

    with open("runmd.sh", 'w') as f:
        f.write('''# Shuoxue Li <sli7@caltech.edu>
# Download sample file from github
filepath="https://raw.githubusercontent.com/LiShuoxue/Ch21C-2023SP/main/Gromacs/Argon-gas/"
for filename in "argon.top" 
do
  curl -O ${filepath}${filename}
done

# change temperature to the arg $1
curl -O ${filepath}md.mdp
cat md.mdp | sed "s/MYTEMP/$1/g" > mdwork.mdp

# Run Molecular dynamics
gmx grompp -f mdwork.mdp -c $2.pdb -p argon.top
gmx mdrun -s topol.tpr -v -c argon_1ns.gro -nice 0

# Post-processing
echo 4  | gmx energy -o Potential.xvg
echo 5  | gmx energy -o Kinetic.xvg
echo 6  | gmx energy -o Total_E.xvg
echo 8  | gmx energy -o Temperature.xvg
echo 10 | gmx energy -o Pressure.xvg
''')

    cnt = 0

    df = pd.DataFrame()
    for tag in ['Potential', 'Kinetic', 'Total_E', 'Temperature', 'Pressure', "Volume"]:
        df[tag] = np.zeros(len(tmps) * len(Vs) * len(Ns))
    for tag in ['Potential', 'Kinetic', 'Total_E', 'Pressure']:
        df[tag+"-std"] = np.zeros(len(tmps) * len(Vs) * len(Ns))
    for tmp in tmps:
        for V in Vs:
            for N in Ns:

                folder = "V-{:.2f}-T-{:.2f}K".format(V, tmp)

                try: os.makedirs(folder)
                except FileExistsError: pass
                os.chdir(folder)
                
                generate_Argon_pdbfile(V, N, "V-{:.2f}-N-{:.0f}".format(V, N))

                df['Volume'][cnt] = V
                df['Temperature'][cnt] = tmp

                subprocess.call(["sh", "../runmd.sh", "{:.2f}".format(tmp), "V-{:.2f}-N-{:.0f}".format(V, N)])
                for tag in ['Potential', 'Kinetic', 'Total_E', 'Pressure']:
                    step, quant = np.loadtxt("{}.xvg".format(tag), comments=['#', '@']).T
                    quant_mean = np.mean(quant[len(step)//3:])
                    quant_std  = np.std(quant[len(step)//3:])

                    print("T = {}, V = {} : {} = ({} ± {}) {}".format(
                        tmp, V, tag, quant_mean, quant_std, units[tag]
                    ))

                    df[tag][cnt] = quant_mean
                    df[tag+'-std'][cnt] = quant_std

                cnt += 1

                os.chdir("..")
                print("-" * 50)

    df.to_csv("results.csv")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-T", "--temperature", help="Select Temperature Range (Start(K), End(K), Step)", nargs='*')
    parser.add_argument("-V", "--volume", help="Select Volume Range (Start(Ang^3), End(Ang^3), Step)", nargs='*')
    parser.add_argument("-N", "--number", help="Select particle numbers (Start, End, Step)")

    args = parser.parse_args()._get_kwargs()

    T_START, T_END, T_STEP = 100, 400, 4
    V_START, V_END, V_STEP = 500000, 500000, 1
    N_START, N_END, N_STEP = 100, 100, 1

    for name, arg in args:
        if name == "temperature":
            if arg is not None:
                T_START, T_END, T_STEP = eval(arg[0]), eval(arg[1]), eval(arg[2])
        if name == "volume":
            if arg is not None:
                V_START, V_END, V_STEP = eval(arg[0]), eval(arg[1]), eval(arg[2])
        if name == "number":
            if arg is not None:
                N_START, N_END, N_STEP = eval(arg[0]), eval(arg[1]), eval(arg[2])

    tmps = np.linspace(T_START, T_END, T_STEP)
    Vs = np.linspace(V_START, V_END, V_STEP)
    Ns = np.linspace(N_START, N_END, N_STEP)

    run_mds(tmps, Vs)
