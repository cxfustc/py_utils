#!/usr/bin/python3

import os
import re
import sys

re_hiden = re.compile (r"^\.")
re_src = re.compile (r"\.c$|\.h$|\.cc$|\.cpp$|\.hh$|\.hpp$|\.py$|\.pl$|\.R$|\.r$|\.f$")

def find_code_file (directory):
	for f in os.listdir(directory):
		if (re_hiden.search(f)):
			continue
		if (os.path.isdir(f)):
			new_dir = directory + "/" + f
			find_code_file (new_dir)
			continue
		if (not re_src.search(f)):
			continue
		cmd = 'vim -c":TOhtml" -c ":wq! ' + out_dir + '/' + directory + '/' + f + '.html" -c ":q" ' + directory + '/' + f
		os.system (cmd)
		cmd = 'wkhtmltopdf ' + out_dir + '/' + directory + '/' + f + '.html ' + out_dir + '/' + directory + '/' + f + '.pdf'
		os.system (cmd)
		cmd = 'rm -f ' + out_dir + '/' + directory + '/' + f + '.html'
		os.system (cmd)

if (len(sys.argv) < 3):
	print ("Usage: python code2pdf <code.dir> <out.dir>")
	sys.exit (-1)

code_dir = os.path.abspath (sys.argv[1])
out_dir = os.path.abspath (sys.argv[2])

if (not os.path.exists(code_dir)):
	print ("'" + code_dir + "' do not exist!")
	sys.exit (-1)

if (not os.path.exists(out_dir)):
	os.makedirs (out_dir)

os.chdir (code_dir)
find_code_file (".")
