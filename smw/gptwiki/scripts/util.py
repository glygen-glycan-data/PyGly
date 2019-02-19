
aamw = {\
            'A': 71.037113848,\
            'C': 103.009185648,\
            'D': 115.026943128,\
            'E': 129.042593208,\
            'F': 147.068414008,\
            'G': 57.021463768,\
            'H': 137.058911944,\
            'I': 113.084064088,\
            'K': 128.094963136,\
            'L': 113.084064088,\
            'M': 131.040485808,\
            'N': 114.042927536,\
            'P': 97.052763928,\
            'Q': 128.058577616,\
            'R': 156.101111152,\
            'S': 87.032028488,\
            'T': 101.047678568,\
            'V': 99.068414008,\
            'W': 186.079313056,\
            'Y': 163.063328648,\
            'n': 0.000000,\
            'c': 18.010560,\
            '+': 1.007825040,\
       }

def peptide_mw(seq):
    mw = aamw['n']+aamw['c']
    for aa in seq:
	mw += aamw[aa]
    return mw

modmw = {\
     '+57.02:C': (+57.021464,''), # Fixed mod
     '+16.00:M': (+15.9994,'Ox') # Var mod
}

def mod_mw(delta,pos,seq):
    if pos == 0:
	aa = '['
	site = aa
    elif pos == (len(seq)+1):
	aa = ']'
	site = aa
    else:
	aa = seq[pos-1]
	site=(aa+str(pos))
    dstr = "%+.2f"%(delta)
    modkey = dstr + ":" + aa
    assert modkey in modmw
    return modmw[modkey][0],site,modmw[modkey][1]

