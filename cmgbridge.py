##Program to fill gaps in time-series files using Joseph's 
# scheme( a spectral method) after filling short gaps with
# a linear fit from neighbors' values.
# 
# [dat,[nogap]]=cmgbridge(dat,nlin,nmaxbr,maxngaps)
# 
# 	dat = data to be bridged. if matrix, bridgeing performed on each column;
#	nogap = 0 if gaps found, 1 if no gaps, optional;
# 	nlin = max gap length to be filled linearly, optional (default =2)
# 	nmaxtbr = max gap length to be filled spectrally, optional (default = 48)
# 	maxngaps = max number of gaps to be filled, optional (default = 1000)
#   outflag - flag whether or not to output to screen. Set to 1 for yes,
#       default is 0 (if left out).
#	
# Adapted from a Fortran program by jpx on 7-15-99
##modified on 12-14-00


import numpy as np

def cmgbridge(dat, nlin, nmaxbr, maxngaps):
    nogap = []
    if dat.ndim == 1:
        dat = dat[:, np.newaxis]
    m, n = dat.shape

    for i in range(n):
        endgap = 0
        notfilled = 0
        thisdata = dat[:, i]
        ngaps, nfirstbad, nlastbad, lgap = cmgidgaps(thisdata, maxngaps)

        if ngaps < 1:
            nogap.append(1)
        else:
            if nlastbad[-1] >= m:
                notfilled += 1
                endgap = 1
            if ngaps > maxngaps:
                notfilled += 1
                ngaps = maxngaps

            redflag = [1]
            for j in range(ngaps):
                if not (endgap and j == ngaps - 1):
                    if lgap[j] <= nlin:
                        nfirstbad[j], nlastbad[j], thisdata = dolinbridge(nfirstbad[j], nlastbad[j], thisdata)

            ngaps2, nfirstbad, nlastbad, lgap = doreseq(ngaps, nfirstbad, nlastbad, lgap)
            if ngaps2 < ngaps:
                redflag.append(0)
            ngaps = ngaps2
            for j in range(ngaps):
                if not (endgap and j == ngaps - 1):
                    nbefore, nafter, notfilled, nflag = dogoodpts(j, m, nmaxbr, ngaps, nfirstbad, nlastbad, lgap, notfilled)
                    if nflag < 1:
                        redflag.append(0)
                        thisdata = dobridgeit(nfirstbad[j], nlastbad[j], lgap[j], nbefore, nafter, thisdata)
            dat[:, i] = thisdata
            nogap.append(all(redflag))

    if len(nogap) > 1:
        return dat, nogap
    else:
        return dat

# function to linearly fill in short gaps in data 
# using the NEIGHBORING good points

def dolinbridge(n1, n2, dat, outflag=0):
    newdat = dat.copy()
    if n1 == 0:
        return -99, -99, newdat
    if n2 == len(dat):
        return -99, -99, newdat
    ng1 = n1 - 1
    ng2 = n2 + 1
    den = ng2 - ng1
    for n in range(n1, n2):
        fact = (n - n1) / den
        newdat[n] = dat[ng1] + fact * (dat[ng2] - dat[ng1])
    return -99, -99, newdat

# function to eliminate linearly-filled gaps from the list of gaps
# START at end of record and work back to first gap that has
# NOT been linearly filled (dolinbridge changed endpoints of
# filled gaps to -99)

def doreseq(ngaps, nfirstbad, nlastbad, lgap):
    indx = np.where(nfirstbad < 0)[0]
    ngaps -= len(indx)
    lgap = np.delete(lgap, indx)
    nfirstbad = np.delete(nfirstbad, indx)
    nlastbad = np.delete(nlastbad, indx)
    return ngaps, nfirstbad, nlastbad, lgap


# Subroutine to look at the length of a data gap and number of 
# good points before and after a data gap to be filled spectrally.  
#	Previously-filled gaps are not good points unless you take a
#	second trip through the program.
#	Subroutine writes message and sets nflag to 1 if it refuses to
#	fill the gap because gap is too long, or to 10 if the lengths of
#	the good points on either end of the gap is too short.
#
#	IMPORTS:
#		igap: gap number
#		npts: number of points in series
#		nmaxbr: max # points to bridge
#		nfirstbad,nlastbad,lgap: arrays of 1st & last points and lengths of gaps
#		ngaps: the number of gaps in the file
# 
#	IMPORTS/EXPORTS
#		notfilled: the # gaps program refuses to fill
#	EXPORTS:
#		nflag: a don't-fill flag: 0=fill 1=gap too long 10=too few good points

def dogoodpts(igap, npts, nmaxbr, ngaps, nfirstbad, nlastbad, lgap, notfilled):
    nflag = 0
    if lgap[igap] > nmaxbr:
        nflag = 1
        notfilled += 1
        nbefore = -999
        nafter = -999
        return nbefore, nafter, notfilled, nflag
    if igap == 0:
        nbefore = nfirstbad[igap] - 1
    else:
        nbefore = nfirstbad[igap] - nlastbad[igap-1] - 1
    if igap == ngaps - 1:
        nafter = npts - nlastbad[igap] - 1
    else:
        nafter = nfirstbad[igap+1] - nlastbad[igap] - 1
    if nbefore < 100 and nbefore < 2 * lgap[igap]:
        nflag = 10
        notfilled += 1
        return nbefore, nafter, notfilled, nflag
    if nafter < 100 and nafter < 2 * lgap[igap]:
        nflag = 10
        notfilled += 1
        return nbefore, nafter, notfilled, nflag
    return nbefore, nafter, notfilled, nflag

def cmgidgaps(dat, maxngaps):
    gaps = np.where(np.isnan(dat))[0]
    ngaps = len(gaps)
    if ngaps > maxngaps:
        ngaps = maxngaps
    nfirstbad = np.zeros(ngaps, dtype=int)
    nlastbad = np.zeros(ngaps, dtype=int)
    lgap = np.zeros(ngaps, dtype=int)
    if ngaps > 0:
        nfirstbad[0] = gaps[0]
        for i in range(1, ngaps):
            nfirstbad[i] = gaps[i]
            nlastbad[i-1] = gaps[i-1]
            lgap[i-1] = gaps[i] - gaps[i-1] - 1
        nlastbad[-1] = gaps[-1]
        lgap[-1] = len(dat) - gaps[-1] - 1
    return ngaps, nfirstbad, nlastbad, lgap

def dobridgeit(n1, n2, lgap, nbefore, nafter, dat):
    newdat = dat.copy()
    if nbefore < 0 or nafter < 0:
        return newdat
    if nbefore < 100 and nbefore < 2 * lgap:
        return newdat
    if nafter < 100 and nafter < 2 * lgap:
        return newdat
    for n in range(n1, n2+1):
        fact = (n - n1) / lgap
        newdat[n] = dat[n1-1-nbefore] + fact * (dat[n2+1+nafter] - dat[n1-1-nbefore])
    return newdat


# Identical to Joseph's UVTBR2
# Modified from a FORTRAN program by jpx on 7-15-99

def dobridge(x, istart, m, mm, lx, ii):
    lxtnd = abs(lx)
    g, p = peflt(x, istart, m, mm)

    if ii == 1:
        g = np.flip(g)
        jstart = istart + m
        jend = jstart + lxtnd - 1
        k = jstart - mm
        for j in range(jstart, jend + 1):
            dotprd = np.sum(x[k:k+mm] * g[:mm])
            x[j] = dotprd
            k += 1
    elif ii == 2:
        jstart = istart - 1
        jend = jstart - lxtnd + 1
        iflag = 0
        k = 0
        frac = 0
        jsave = 0
        for j in range(jstart, jend - 1, -1):
            k += 1
            dotprd = np.sum(x[j+1:j+mm+1] * g[:mm])
            if not np.isnan(x[j]):
                if iflag == 0:
                    iflag = 1
                    jsave = j
                    if lx < 0:
                        dfrac = 0.333333333
                    else:
                        if lxtnd == 1:
                            dfrac = 0.5
                        else:
                            dfrac = 1 / (lxtnd - k + 1)
                else:
                    frac += dfrac
                scrtch = frac * x[j] + (1 - frac) * dotprd
            x[j] = dotprd
        if jsave != 0:
            k = jstart - jsave
            for j in range(jsave, jend - 1, -1):
                k += 1
                x[j] = scrtch[k]

    return x


# 	BOTTERO'S SUBROUTINE BASED ON ANDERSON'S EQUATIONS IN
# 	VOL 39 NO 1 OF Geophysics 1974, p69-72.
#   Modified from a FORTRAN program by jpx on 7-15-99
def peflt(x, istart, n, mmax):
    iend = istart + n - 1
    sum1 = np.sum(x[istart:iend+1] ** 2)
    p = sum1 / n

    m = 1
    nmm = n - m
    b = x[istart:istart+nmm]
    bp = x[istart+1:istart+nmm+1]

    while m <= mmax:
        if m > 1:
            nmm = n - m
            for i in range(nmm):
                b[i] = b[i] - aold[m-1] * bp[i]
                bp[i] = bp[i+1] - aold[m-1] * b[i+1]

        sum1 = np.sum(b[:nmm] * bp[:nmm])
        sum2 = np.sum(b[:nmm]**2 + bp[:nmm]**2)
        if sum2 != 0:
            a = 2 * sum1 / sum2
        else:
            a = 0
        p *= (1 - a**2)
        if m > 1:
            mm1 = m - 1
            for k in range(mm1):
                mmk = m - k - 1
                a[k] = aold[k] - a * aold[mmk]
        aold = a.copy()
        m += 1

    g = a
    return g, p

def dobridgeit(nfirstbad, nlastbad, lgap, nbefore, nafter, dat):
    lxtnd = round(0.75 * lgap)
    if lxtnd <= 0:
        lxtnd = 1
    elif lgap == 2:
        lxtnd = -2

    nbefore = min(nbefore, 200)
    istart = nfirstbad - nbefore
    nb2 = int(nbefore // 2)
    ii = 1
    dat = dobridge(dat, istart, nbefore, nb2, lxtnd, ii)

    istart = nlastbad + 1
    nafter = min(nafter, 200)
    na2 = int(nafter // 2)
    ii = 2
    dat = dobridge(dat, istart, nafter, na2, lxtnd, ii)

    return dat

