import numpy as np
from scipy.signal import csd, welch

def reynolds_stress(au, av, aw, bu, bv, bw, fs, bn, jd):
    DirA = 0
    ialign = 0
    zref = 0.35

    xd = au - bu
    yd = av - bv
    zd = aw - bw

    cxw = np.cov(np.column_stack((xd, zd)), ddof=1)
    rsu = -0.5 * cxw[0, 1]
    cyw = np.cov(np.column_stack((yd, zd)), ddof=1)
    rsv = -0.5 * cyw[0, 1]

    rss = {
        'rsu': rsu,
        'rsv': rsv
    }

    rs_mag, rs_dir = pcoord(rsu, rsv)
    rss['mag'] = rs_mag
    rss['dir'] = rs_dir

    sam, dam = pcoord(np.mean(au), np.mean(av))
    sbm, dbm = pcoord(np.mean(bu), np.mean(bv))
    sabm, dabm = pcoord(np.mean(np.concatenate((au, bu))), np.mean(np.concatenate((av, bv))))

    sa, da = pcoord(au, av)
    sb, db = pcoord(bu, bv)

    da_r = da - dabm + 90
    db_r = db - dabm + 90
    u_ar, v_ar = xycoord(sa, da_r)
    u_br, v_br = xycoord(sb, db_r)

    sarm, darm = pcoord(np.mean(u_ar), np.mean(v_ar))
    sbrm, dbrm = pcoord(np.mean(u_br), np.mean(v_br))

    xd = u_ar - u_br
    yd = v_ar - v_br
    zd = aw - bw

    cxw = np.cov(np.column_stack((xd, zd)), ddof=1)
    rsur = -0.5 * cxw[0, 1]
    cyw = np.cov(np.column_stack((yd, zd)), ddof=1)
    rsvr = -0.5 * cyw[0, 1]

    dof = sabm * len(xd) / (fs * zref)
    e1x = 0.5 * np.sqrt(np.var(xd) * np.var(zd)) / np.sqrt(dof)
    e1y = 0.5 * np.sqrt(np.var(yd) * np.var(zd)) / np.sqrt(dof)
    e1 = np.sqrt(e1x ** 2 + e1y ** 2)

    e2x = rsur * 0.5 * np.sqrt((1 + (cxw[0, 0] * cxw[1, 1]) / cxw[0, 1] ** 2) / dof)
    e2y = rsvr * 0.5 * np.sqrt((1 + (cyw[0, 0] * cyw[1, 1]) / cyw[0, 1] ** 2) / dof)
    e2 = np.sqrt(e2x ** 2 + e2y ** 2)

    rss['e2x'] = e2x
    rss['e2y'] = e2y
    rss['e2'] = e2
    rss['dof'] = dof

    pct = -1
    n = len(xd)
    n4 = int(n / 4)
    xxz, lags = np.correlate(xd, zd, 'full')
    xyz, lags = np.correlate(yd, zd, 'full')
    ax, lags = np.correlate(xd, 'full')
    ay, lags = np.correlate(yd, 'full')
    az, lags = np.correlate(zd, 'full')
    lag0 = np.where(lags == 0)[0][0]

    timescale = zref * fs / sabm
    print(f'time scale (lags) = {timescale}')

    nseg = 16
    nfft = len(xd) // nseg
    Pxz, ff = welch(xd, zd, fs, nperseg=nfft, window='hanning', return_onesided=False)
    Pyz, ff = welch(yd, zd, fs, nperseg=nfft, window='hanning', return_onesided=False)

    nseg = 1
    nfft = len(xd) // nseg
    Pxz, ff = welch(xd, zd, fs, nperseg=nfft, window='boxcar', return_onesided=False)
    Pyz, ff = welch(yd, zd, fs, nperseg=nfft, window='boxcar', return_onesided=False)
    dff = ff[1] - ff[0]
    dffp = dff

    rss['pct'] = pct

    return rss

def pcoord(u, v):
    mag = np.sqrt(u ** 2 + v ** 2)
    dir = np.arctan2(v, u) * 180 / np.pi
    return mag, dir

def xycoord(r, theta):
    x = r * np.cos(theta * np.pi / 180)
    y = r * np.sin(theta * np.pi / 180)
    return x, y

