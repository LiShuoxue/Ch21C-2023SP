import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-n", "--nsite", help="Number of site")
parser.add_argument("-U", "--interaction", help="Hubbard U/t")
parser.add_argument("-d", "--bond-dimension", help="DMRG bond dimension")
args = parser.parse_args()._get_kwargs()
for name, arg in args:
    print(name)
    print(arg)
    print(type(arg))