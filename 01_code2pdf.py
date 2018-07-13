#!/usr/bin/python3

import os
import re
import sys

re_hiden = re.compile (r"^\.")
re_src = re.compile (r"\.c$|\.h$|\.cc$|\.cpp$|\.hh$|\.hpp$|\.py$|\.pl$|\.R$|\.r$|\.f$|\.vim$|\.css$")

def find_code_file (directory):
	for f in os.listdir(directory):
		if (re_hiden.search(f)):
			continue
		path = directory + "/" + f
		if (os.path.isdir(path)):
			find_code_file (path)
			continue
		if (not re_src.search(f)):
			continue
		cur_out_dir = out_dir + '/' + directory
		if (not os.path.exists(cur_out_dir)):
			os.makedirs (cur_out_dir)

		cmd = 'vim -c":TOhtml" -c ":wq! ' + cur_out_dir + '/' + f + '.html" -c ":q" ' + directory + '/' + f
		os.system (cmd)
		cmd = 'wkhtmltopdf ' + cur_out_dir + '/' + f + '.html ' + cur_out_dir + '/' + f + '.pdf'
		os.system (cmd)
		cmd = 'rm -f ' + cur_out_dir + '/' + f + '.html'
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
