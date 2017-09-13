#!/usr/bin/env python

import sys, getopt
import math as ma        

def calcNvecs(box_length, cut_off, precision):
    TWOpi = 2.0 * ma.pi
    alpha = ma.sqrt(precision) / cut_off
    h_ewald = 2.0 * precision / cut_off
    basis_length = box_length
    k_max = int(h_ewald * basis_length / TWOpi) + 1
    hcut_sq = h_ewald * h_ewald

    nvecs = 0
    for nz in range(-k_max,k_max+1):
        for ny in range(-k_max,k_max+1):
            for nx in range(k_max+1):
                if (nx == 0 and ny == 0 and nz == 0):
                    continue

                hx = TWOpi * nx / basis_length 
                hy = TWOpi * ny / basis_length 
                hz = TWOpi * nz / basis_length 
                hsq = hx*hx + hy*hy + hz*hz
                if (hsq < hcut_sq):
                    nvecs += 1
    return nvecs
def main(argv=None):
    if argv is None:
        argv = sys.argv

    try:
        opts, args = getopt.getopt(argv[1:], "hl:c:t:v:", 
                     ["help", "box length=", "cut_off=", "tolerance=", "nvecs="])
    except getopt.error, msg:
        print msg
        print "for help use --help"
        return 2

    flags = []
    n_config = 1
    calc_tot = False
    for opt, arg in opts:
        if opt == '-h':
            print "calcNvecs.py -l <box_length> -c <cut_off> -t <tolerance> -v <nvecs>"
            return 1
        elif opt == '-l':
            box_length = float(arg)
        elif opt == '-c':
            cut_off = float(arg)
        elif opt == '-t':
            precision = -ma.log(float(arg))
        elif opt == '-v':
            nvecs = int(arg)
            calc_tot = True

    if (calc_tot):
        precision = 2.3
        t_nvecs = 0
        while (t_nvecs != nvecs):
            if (t_nvecs <= nvecs):
                precision *= 1.1
            else:
                precision *= 0.1

            t_nvecs = calcNvecs(box_length, cut_off, precision)

        print 'tolerance= ', ma.exp(-precision)

    else:
        print 'nvecs= ', calcNvecs(box_length, cut_off, precision)
    
if __name__ == '__main__':
    sys.exit(main())

# "The greatest happiness is to know the source of unhappiness." -Fyodor Dostoevsky
