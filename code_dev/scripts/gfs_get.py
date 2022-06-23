# gfs_get.py
## Eden McEwen
## 10.2021
## Takes in a RUN file and pulls the gfs data for those dates

from pipeline.code.file_reader import *
from urllib.request import urlretrieve

prefix = "http://mkwc.ifa.hawaii.edu/forecast/mko/models/gfs/"
wspd_suffix = ".wspd.txt"
wdir_suffix = ".wdir.txt"

def do_url(date, out_d):
    entry = ["00", "06", "12", "18"]
    for e in entry:
        file = f"{prefix}gfs.{date}{e}{wspd_suffix}"
        try:
            urlretrieve(f"{prefix}{date}{e}/gfs.{date}{e}{wspd_suffix}", 
                           filename=f"{out_d}gfs.{date}{e}{wspd_suffix}")
            urlretrieve(f"{prefix}{date}{e}/gfs.{date}{e}{wdir_suffix}", 
                           filename=f"{out_d}gfs.{date}{e}{wdir_suffix}")
        except:
            print("Not Found: ", file)


if __name__ == '__main__':
    """
    arg1: RUN file 
    arg2: out_dir
    """
    run_file = sys.argv[1]
    out_dir = sys.argv[2]
    
    entries = read_file(run_file)
    
    for d in entries:
        print(d)
        do_url(d, out_dir)


