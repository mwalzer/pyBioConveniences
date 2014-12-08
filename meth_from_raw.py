#!/usr/bin/env python
"""
	Read method file name from raw file - not recursive
"""
__author__ = 'walzer'

import sys
import argparse
import re
from os import listdir
from os.path import isfile, isdir, join
import pickle

def access_meth(rf):
	if isfile(rf):
		met_regex = None
		with open(rf, 'rb') as f:
			met_regex = re.findall( b'\x00\x43\x00\x3A\x00(.*?)\x00\x6D\x00\x65\x00\x74\x00\x68\x00\x00\x00\x00', f.read() )
		try:
			method_name = met_regex[0].replace('\x00','').split('\\')[-1] + 'meth'
			return method_name
		except:
			pass
	return None


def main():
	"""
	manage input & output, call ascii access function
	"""
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument('-i','--input', help='Raw file or directory with files to read raw from.', required=True)
	parser.add_argument('-p','--pickle', help='Pickle for result dictionary, no console output.', required=False)
	print parser.description
	if len(sys.argv)<=1:
		parser.print_help()
		sys.exit(1)
	args = vars(parser.parse_args())
	if isfile(args['input']) and (args['input'].endswith('.RAW') or args['input'].endswith('.raw') ):
		m = access_meth(args['input'])
		if args['pickle']:
			pickle.dump( {args['input']:m}, open( args['pickle'], "wb" ) )
		else:
			print args['input'], ":", m
	elif isdir(args['input']):
		onlyrawfiles = [ join(args['input'],f) for f in listdir(args['input']) if ( isfile(join(args['input'],f)) and (f.endswith('.RAW') or f.endswith('.raw') ) and (not f.startswith('~')))]
		m = {f: access_meth(f) for f in onlyrawfiles}
		if args['pickle']:
			pickle.dump( m, open( args['pickle'], "wb" ) )
		else:
			for p in m:
				print p, ":", m[p]
	else:
		parser.print_help()

if __name__ == "__main__":
    main()