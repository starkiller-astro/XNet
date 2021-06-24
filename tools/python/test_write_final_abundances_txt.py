from write_final_abundances_txt import write_final_abundances
import sys
import argparse

parser = argparse.ArgumentParser(description="Takes final abuncances from a set of ts files and combines in one output file.")
parser.add_argument('path', type=str,  nargs=1, help="path to ts output files, e.g.: /results/ts* ")
parser.add_argument('-o','--output', type=str,  nargs=1, help="output filename",default='final_abundances.txt')
parser.add_argument('-a', '--append', action='store_const', const=True, help="flag to append to existing file.",default=False)

args = parser.parse_args()
path=args.path[0]
output = args.output[0]
append = args.append

write_final_abundances(path, file_name = output, appending = append, print_positions=False)

