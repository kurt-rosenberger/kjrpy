import numpy as np
from scipy.integrate import quadrature

def dissipation_rate(uvw, fs, nfft=512, iplot=0):
    # Units in output assume units in input are MKS
    alpha = 1.5  # Kolomogorov constant

    ss = uvwstats(uvw[:, 0], uvw[:, 1])
    S = ss.S
    sd = ss.sd1
    phir = ss.phir

    noverlap = nfft // 2
    window = np.hamming(nfft)
    nf = (nfft // 2) + 1
    U, f = welch(uvw[:, 0], window, noverlap, nfft, fs)
    U = U[1:nf, 0]
    V, f = welch(uvw[:, 1], window, noverlap, nfft, fs)
    V = V[1:nf, 0]
    W, f = welch(uvw[:, 2], window, noverlap, nfft, fs)
    W = W[1:nf, 0]
    f = f[1:nf, 0]

    # dissipation rate estimates
    twopi = 2 * np.pi
    IA = bhi(sd / S, phir)
    # Frequency range that is used for inertial-range fit
    fm = np.where((f >= 0.6) & (f <= 2.5))[0]
    fhn = np.where(f >= 3.5)[0]  # assume noise at f > 3.5 Hz
    uvnoise = np.mean(U[fhn, 0] + V[fhn, 0])
    N = uvnoise * np.ones_like(U)
    Puv = (twopi * f) ** (5 / 3) * np.maximum((U + V - N) / twopi, np.eps * np.ones_like(U))
    Pww = (twopi * f) ** (5 / 3) * W / twopi
    euv = np.mean(Puv[fm]) * (55 / (21 * alpha)) * (S ** (-2 / 3)) * (1 / IA)
    eww = np.mean(Pww[fm]) * (55 / (12 * alpha)) * (S ** (-2 / 3)) * (1 / IA)
    euv_sd = np.std(Puv[fm]) * (55 / (21 * alpha)) * (S ** (-2 / 3)) * (1 / IA)
    eww_sd = np.std(Pww[fm]) * (55 / (12 * alpha)) * (S ** (-2 / 3)) * (1 / IA)

    return {'euv': euv, 'eww': eww, 'euv_sd': euv_sd, 'eww_sd': eww_sd, 'uvnoise': uvnoise, 'IA': IA, 'ss': ss}

def bhi(sv, phi):
    F = f'((x**2 - 2*({1/sv})*np.cos({phi})*x + ({1/sv**2}))**(1/3)*np.exp(-0.5*x**2))'
    f = lambda x: eval(F)
    I, _ = quadrature(f, -5, 5)
    return (1 / np.sqrt(2 * np.pi)) * (sv ** (2 / 3)) * I

def uvwstats(u, v, sf=0, ipost=0):
    mu = np.mean(u)
    mv = np.mean(v)
    C = np.cov(u, v)
    V, D = np.linalg.eig(C)

    x1 = [0.5 * np.sqrt(D[0, 0]) * V[0, 0], -0.5 * np.sqrt(D[0, 0]) * V[0, 0]]
    y1 = [0.5 * np.sqrt(D[0, 0]) * V[1, 0], -0.5 * np.sqrt(D[0, 0]) * V[1, 0]]
    x2 = [0.5 * np.sqrt(D[1, 1]) * V[0, 1], -0.5 * np.sqrt(D[1, 1]) * V[0, 1]]
    y2 = [0.5 * np.sqrt(D[1, 1]) * V[1, 1], -0.5 * np.sqrt(D[1, 1]) * V[1, 1]]
    mspd, mdir = pcoord(mu, mv)
    l1, az1 = pcoord(x1[0], y1[0])
    l2, az2 = pcoord(x2[0], y2[0])
    if l1 < l2:
        l1, az1, l2, az2 = l2, az2, l1, az1
    sd1 = 2 * l1
    sd2 = 2 * l2

    return Bunch(ubar=mu, vbar=mv, S=mspd, az0=mdir, sd1=sd1, sd2=sd2, az1=az1, az2=az2, phid=az1 - mdir, phir=np.abs(az1 - mdir) * np.pi / 180)

def pcoord(x, y):
    r = np.sqrt(x ** 2 + y ** 2)
    q1 = (x > 0) & (y >= 0)
    q2 = (x < 0) & (y >= 0)
    q3 = (x < 0) & (y < 0)
    q4 = (x > 0) & (y < 0)
    north = (x == 0) & (y > 0)
    south = (x == 0) & (y < 0)
    zip_ = (x == 0) & (y == 0)
    small = (x == 0)
    xs = x.copy()
    xs[small] = np.ones_like(small) * np.finfo(float).eps
    rad = 180 / np.pi
    ang = rad * np.arctan(y / xs)
    az = q1 * (90 - ang) + q2 * (270 - ang) + q3 * (270 - ang) + q4 * (90 - ang)
    if len(north) > 0:
        az[north] = np.ones_like(north) * 360
    if len(south) > 0:
        az[south] = np.ones_like(south) * 180
    if len(zip_) > 0:
        az[zip_] = np.zeros_like(zip_)
    if np.any(az >= 360):
        i360 = np.where(az >= 360)[0]
        az[i360] = az[i360] - 360
    return r, az

def welch(x, window, noverlap, nfft, fs):
    f, Pxx = signal.welch(x, fs, window=window, nperseg=nfft, noverlap=noverlap)
    return Pxx, f

class Bunch:
    def __init__(self, **kwds):
        self.__dict__.update(kwds)

