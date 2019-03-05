#!/usr/bin/env python3

import sys, argparse, os
sys.path.append('../../../trunk/tools/python')

import find_difference as fd

parser = argparse.ArgumentParser(description='Compare two binary output files from XNet')
parser.add_argument('file1', type=argparse.FileType('r'), help='Name of 1st binary file')
parser.add_argument('file2', type=argparse.FileType('r'), help='Name of 2nd binary file')
args = parser.parse_args()

fd.compare_final(args.file1.name,args.file2.name)