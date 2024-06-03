import numpy as np

def cmgdataclean(gin, threshold=1e34):
    indx = np.where(gin >= threshold)[0]
    gin[indx] = np.nan
    indx2 = np.where(gin <= threshold)[0]
    gin[indx2] = np.nan
    gout = gin
    return gout

def cmgpca(u, v, onoff):
    indx = np.where(np.isnan(u) == 0)[0]
    u = u[indx]
    v = v[indx]
    indx = np.where(np.isnan(v) == 0)[0]
    u = u[indx]
    v = v[indx]
    mu = np.mean(u)
    mv = np.mean(v)
    covar = np.cov(u, v)
    eigvec, eigval = np.linalg.eig(covar)
    spd, direc = cmguv2spd(mu, mv)

    if eigval[1, 1] > eigval[0, 0]:
        cols = [1, 0]
    else:
        cols = [0, 1]
    majr = np.sqrt(eigval[cols[0], cols[0]])
    junk, majaz = cmguv2spd(eigvec[0, cols[0]], eigvec[1, cols[0]])
    minr = np.sqrt(eigval[cols[1], cols[1]])
    junk, minaz = cmguv2spd(eigvec[0, cols[1]], eigvec[1, cols[1]])
    x1, y1 = cmgspd2uv(majr, majaz)
    x2, y2 = cmgspd2uv(minr, minaz)

def cmgpolar(u, v):
    u1 = u.copy()
    v1 = v.copy()
    u1 = u1[np.where(np.isnan(u1) == False)]
    v1 = v1[np.where(np.isnan(v1) == False)]
    if np.all(u1 >= 0) and np.all(v1 >= 0):
        u = u
    else:
        u, v = cmguv2spd(u, v)
    v = v * np.pi / 180
    umax = np.max(u)
    umin = np.min(u)
    udiff = umax - umin
    if udiff < 1000:
        umax2 = 100 * np.ceil(umax / 100)
    if udiff < 100:
        umax2 = 10 * np.ceil(umax / 10)
    if udiff < 10:
        umax2 = np.ceil(umax)
    import matplotlib.pyplot as plt
    plt.polar([v, 0], [u, umax2], '+')
    plt.view(90, -90)

def cmgrotate(east, north, theta):
    newx = east * np.cos(theta) + north * np.sin(theta)
    newy = -east * np.sin(theta) + north * np.cos(theta)
    return newx, newy

def cmgspd2uv(spd, direc):
    spd = cmgdataclean(spd)
    direc = cmgdataclean(direc)
    north, east = np.pol2cart(direc * np.pi / 180, spd)
    return north, east

def cmguv2spd(east, north):
    east = cmgdataclean(east)
    north = cmgdataclean(north)
    direc, spd = np.cart2pol(north, east)
    direc = direc * 180 / np.pi
    indx = np.where(direc < 0)[0]
    direc[indx] = direc[indx] + 360
    return spd, direc

