import numpy as np

def cmgeof(x, *args):
    opslen = len(args) % 2
    m, n = x.shape
    num = n
    donorm = 0
    covyes = 0

    while args:
        s1 = args[0]
        if s1[:4].lower() == 'mode':
            num = args[1]
        elif s1[:4].lower() == 'norm':
            donorm = args[1]
        elif s1[:4].lower() == 'covy':
            covyes = args[1]
        args = args[2:]

    if num > n or num < 1:
        num = n

    if covyes:
        if not np.array_equal(np.tril(x), np.rot90(np.flipud(np.triu(x)), -1)) and np.isreal(x):
            raise ValueError('The input covariance matrix must be symmetric!')
        z = x
        c = z
    else:
        if donorm:
            if np.isreal(x):
                z = signal.detrend(x, axis=0)
                z = (z - np.mean(z, axis=0)) / np.std(z, axis=0)
            else:
                xr = np.real(x)
                xi = np.imag(x)
                zr = signal.detrend(xr, axis=0)
                zi = signal.detrend(xi, axis=0)
                zr = (zr - np.mean(zr, axis=0)) / np.std(zr, axis=0)
                zi = (zi - np.mean(zi, axis=0)) / np.std(zi, axis=0)
                z = zr + 1j * zi
        else:
            z = signal.detrend(x, axis=0)
        c = np.cov(z.T)

    vv = np.diag(c)
    totalv = np.sum(vv)

    eigvec, eigval = np.linalg.eig(c)

    mv = np.diag(eigval)
    jj = np.argsort(mv)
    mv = np.flip(mv[jj])
    jj = np.flip(jj)
    eigval = mv

    eigvec = eigvec[:, jj]

    sig = np.sign(eigvec)
    ampsphase = np.full_like(sig, np.nan)

    if covyes:
        ampt = np.full((1, len(mv)), np.nan)
        amptr = ampt
        ampti = ampt
        print('\nTemporal amplitude can not be computed when the input is a convariance matrix.\n')
    else:
        ampt = z @ eigvec
        amptr = np.real(ampt) / np.sqrt(np.ones((m, 1)) * mv)
        ampti = np.imag(ampt) / np.sqrt(np.ones((m, 1)) * mv)

    mv = mv * np.ones((1, len(mv)))
    sv = mv * (eigvec * np.conj(eigvec))
    amps = np.sqrt(2 * sv)
    vv = vv * np.ones((1, len(vv)))
    sv = 100 * sv / vv
    mv = 100 * np.diag(mv) / totalv

    if num < n:
        ampt = ampt[:, :num]
        amps = amps[:, :num]
        sv = sv[:, :num]
        mv = mv[:num]
        eigvec = eigvec[:, :num]
        sig = sig[:, :num]
        ampsphase = ampsphase[:, :num]

    if np.isreal(x):
        amps = amps * sig
    else:
        ampsphase = 180 / np.pi * np.angle(eigvec)

    return ampt, amps, sv, mv, eigvec, sig, ampsphase

