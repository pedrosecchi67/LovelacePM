import numpy as np

def read_afl(afl, ext_append=False, header_lines=1): #read arfoil data points from Selig format file
    if ext_append:
        infile=open(afl+'.dat', 'r')
    else:
        infile=open(afl, 'r')
    aflpts=[]
    alltext=infile.read()
    lines=alltext.split('\n')
    for i in range(header_lines, len(lines)):
        linelist=lines[i].split()
        aflpts+=[[float(linelist[0]), float(linelist[1])]]
    return np.array(aflpts)