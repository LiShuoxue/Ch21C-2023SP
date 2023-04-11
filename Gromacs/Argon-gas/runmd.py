import os, subprocess, argparse
import pandas as pd
import numpy as np

units = {
    "Potential": "kJ/mol",
    "Kinetic": "kJ/mol",
    "Total_E": "kJ/mol",
    "Temperature": "K",
    "Pressure": "Bar",
    "Volume": "L/mol"
}



def generate_Argon_pdbfile(V, N, filename=None):

    a = V**(1/3)

    Nsup = int(N**(1/3)) + 1

    if filename is None: filename = "argon-start-N{}-V{}.pdb".format(N, V)

    with open(filename+".pdb", "w") as f:
        f.write("HEADER    Argon\n")
        f.write("REMARK    THIS IS A SIMULATION BOX\n")
        f.write("CRYST1 {:.6f} {:.6f} {:.6f}  90.00  90.00  90.00 P 1  1\n".format(a,a,a))
        f.write("MODEL        1\n")

        idx = 0

        while idx < 600:
            for i in range(Nsup):
                for j in range(Nsup):
                    for k in range(Nsup):

                        idx += 1
                        f.write("ATOM{:>7}  Ar   Ar{:>6}{:>12.3f}{:>8.3f}{:>8.3f}  1.00  0.00\n".format(idx+1, idx+1, i*a/Nsup, j*a/Nsup, k*a/Nsup))
        
        """
        for i in range(N):
            f.write("ATOM{:>7}  Ar   Ar{:>6}{:>12.3f}{:>8.3f}{:>8.3f}  1.00  0.00\n".format(i+1, i+1, np.random.random()*a, np.random.random()*a, np.random.random()*a))
        """
        f.write("TER\n")
        f.write("ENDMDL\n")
        f.close()

def run_mds(tmps, Vs, Ns):

    with open("runmd.sh", 'w') as f:
        f.write('''# Shuoxue Li <sli7@caltech.edu>
# Download sample file from github
filepath="https://raw.githubusercontent.com/LiShuoxue/Ch21C-2023SP/main/Gromacs/Argon-gas/"
for filename in "argon-template.top" 
do
  curl -O ${filepath}${filename}
done

cat argon-template.top | sed "s/NUMBER/$3/g" > argon.top

# change temperature to the arg $1
curl -O ${filepath}md.mdp
cat md.mdp | sed "s/MYTEMP/$1/g" > mdwork.mdp

# Run Molecular dynamics
gmx grompp -f mdwork.mdp -c V-$2-N-$3.pdb -p argon.top
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
    for tag in ['Potential', 'Kinetic', 'Total_E', 'Temperature', 'Pressure', "Volume", "Number"]:
        df[tag] = np.zeros(len(tmps) * len(Vs) * len(Ns))
    for tag in ['Potential', 'Kinetic', 'Total_E', 'Pressure']:
        df[tag+"-std"] = np.zeros(len(tmps) * len(Vs) * len(Ns))
    for tmp in tmps:
        for Vm in Vs:
            for N in Ns:

                V = (Vm * 1E-3) * N / 6.02214076e+23 * 1E30 # Vm (L/mol) -> V (Ang^3)
                N = int(N)

                folder = "V-{:.2f}-T-{:.2f}K-N-{:.0f}".format(Vm, tmp, N)

                try: os.makedirs(folder)
                except FileExistsError: pass
                os.chdir(folder)
                
                generate_Argon_pdbfile(V, N, "V-{:.2f}-N-{:.0f}".format(Vm, N))

                df['Volume'][cnt] = Vm
                df['Temperature'][cnt] = tmp
                df['Number'][cnt] = N

                subprocess.call(["sh", "../runmd.sh", "{:.2f}".format(tmp), "{:.2f}".format(Vm), "{:.0f}".format(N)])
                for tag in ['Potential', 'Kinetic', 'Total_E', 'Pressure']:
                    step, quant = np.loadtxt("{}.xvg".format(tag), comments=['#', '@']).T
                    quant_mean = np.mean(quant[len(step)//3:])
                    quant_std  = np.std(quant[len(step)//3:])

                    print("T = {} K, V = {} L/mol, N = {}: {} = ({} ± {}) {}".format(
                        tmp, Vm, N, tag, quant_mean, quant_std, units[tag]
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
    parser.add_argument("-V", "--volume", help="Select Volume per mol Range (Start(L/mol), End(L/mol), Step)", nargs='*')
    parser.add_argument("-N", "--number", help="Select particle numbers (Start, End, Step)",  nargs='*')

    args = parser.parse_args()._get_kwargs()

    T_START, T_END, T_STEP = 100, 400, 4
    V_START, V_END, V_STEP = 0.2, 1, 5
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
                N_START, N_END, N_STEP = int(eval(arg[0])), int(eval(arg[1])), int(eval(arg[2]))

    tmps = np.linspace(T_START, T_END, T_STEP)
    Vs = np.linspace(V_START, V_END, V_STEP)
    Ns = np.linspace(N_START, N_END, N_STEP)

    run_mds(tmps, Vs, Ns)
