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

def run_mds(tmps):
    
    logger = []

    df = pd.DataFrame()
    for tag in ['Potential', 'Kinetic', 'Total_E', 'Temperature', 'Pressure']:
        df[tag] = np.zeros(len(tmps))
    for i, tmp in enumerate(tmps):

        folder = "T-{:.2f}K".format(tmp)

        try: os.makedirs(folder)
        except FileExistsError: pass
        os.chdir(folder)

        subprocess.run(["sh ../runmd.sh", "{:.2f}".format(tmp)])
        for tag in ['Potential', 'Kinetic', 'Total_E', 'Temperature', 'Pressure']:
            step, quant = np.loadtxt("{}.xvg".format(tag), comments=['#', '@']).T
            quant_mean = np.mean(quant)
            quant_std  = np.std(quant)

            logger.append("Point {} : {} = ({} ± {}) {}".format(
                i+1, tag, quant_mean, quant_std, units[tag]
            ))

            df[tag][i] = quant_mean

        os.chdir("..")
        logger.append("----------------------------------------")

    for strs in logger: print(strs)
    df.to_csv("results.csv")

