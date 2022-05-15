import urllib.request
import re
import os

html = urllib.request.urlopen(
    "https://github.com/TLCFEM/suanPan/releases").read()

version = re.search(r'suanPan-v[0-9]\.[0-9]\.?[0-9]?', html.decode('utf-8'))

fn = './.latest.version.sp'

if os.path.exists(fn):
    os.remove(fn)

with open(fn, 'w') as f:
    f.write(version.group())

f.close()
