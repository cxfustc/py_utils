#!/usr/bin/python

import re
import sys
import gzip

re_info = re.compile (r"^##")
re_head = re.compile (r"^#")

if (len(sys.argv) < 2):
  print ("Usage: python vcf_ana <vcf>")
  sys.exit (-1)

vcf = sys.argv[1]

with gzip.open (vcf) as f:
  lines = f.read ();

for line in lines:
  if (re_info.search(line)):
    continue
  if (re_head.search(line)):
    print (line)
    break;
