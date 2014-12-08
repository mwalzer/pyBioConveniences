__author__ = 'walzer'

#! /usr/bin/python

import sys
import argparse
import re
from os import listdir
from os.path import isfile, isdirectory, join
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


def __main__():
	parser = argparse.ArgumentParser(description='Read method file name from raw file - not recursive')
	parser.add_argument('-in','--input', help='Raw file or directory with files to read raw from.', required=True)
	parser.add_argument('-p','--pickle', help='Pickle for result dictionary, no console output.', required=False)
	if len(sys.argv)==1:
		parser.print_help()
		sys.exit(1)
	args = vars(parser.parse_args())
	if not args['in']:
		parser.print_help()
		
	if isfile(args['in']) and (args['in'].endswith('.RAW') or args['in'].endswith('.raw') ):
		m = access_meth(args['in'])
		if args['p']:
			pickle.dump( {args['in']:m}, open( args['p'], "wb" ) )
		else:
			print args['in'], ":", m
	else if isdirectory(args['in']):
		onlyrawfiles = [ f for f in listdir(args['in']) if ( isfile(join(args['in'],f)) and (f.endswith('.RAW') or f.endswith('.raw') ) and (not f.startswith('~')))]
		m = {f: access_meth(f) for f in onlyrawfiles}
		if args['p']:
			pickle.dump( m, open( args['p'], "wb" ) )
		else:
			for p in m:
				print p, ":", m[p]
	else:
		parser.print_help()

