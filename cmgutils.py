import mat2py as mp
from mat2py.core import *


def cmgdataclean(gin, threshold):
    if nargin < 1:
        help(mfilename)
        return gout
    if nargin < 2:
        threshold = 1e34
    indx = find(gin >= threshold)
    gin[I[indx]] = copy(nan)
    indx2 = find(gin <= threshold)
    gin[I[indx2]] = copy(nan)
    gout = copy(gin)
    return gout
    return gout


def cmguv2spd(east, north):
    if nargin < 2:
        help(mfilename)
        return spd, direc
    if any(size(east) - size(north)):
        fprintf("\nTwo input arguments must be the same size.\n")
        return spd, direc
    east = cmgdataclean(east)
    north = cmgdataclean(north)
    direc, spd = cart2pol(north, east)
    direc = mrdivide(direc * 180, pi)
    indx = find(direc < 0)
    direc[I[indx]] = direc(indx) + 360
    return spd, direc
    return spd, direc


def cmgspd2uv(spd, direc):
    if nargin < 2:
        help(mfilename)
        return east, north
    if any(size(spd) - size(direc)):
        fprintf("\nTwo input arguments must be the same size.\n")
        return east, north
    spd = cmgdataclean(spd)
    direc = cmgdataclean(direc)
    north, east = pol2cart((direc @ M[pi]) / 180, spd)
    return east, north
    return east, north


def cmglocalk(pd, dep):
    omega = (2 * pi) / pd
    g = 9.81
    a = mrdivide((omega**2) * dep, g)
    l0 = 1.56 * (pd**2)
    x0 = ((2 * pi) @ M[dep]) / l0
    eps = 1000
    i = 0
    while eps > 1e-8:
        i = i + 1
        fx = (x0 * tanh(x0)) - a
        fxd = (x0 * (1 - (tanh(x0) * tanh(x0)))) + tanh(x0)
        x = x0 - (fx / fxd)
        eps = abs(1 - (x / x0))
        x0 = copy(x)
        if i >= 100:
            break

    k = x / dep
    return k
    return k


def cmgrotate(data, varargin):
    if nargin < 2:
        help(mfilename)
        return varargout
    if length(varargin) > 1:
        east = copy(data)
        north = varargin[I[1]]
        theta = (varargin[I[2]] @ M[pi]) / 180
        newx = (east * cos(theta)) + (north * sin(theta))
        newy = ((-east) * sin(theta)) + (north * cos(theta))
        if _not(isequal(length(east), length(north))):
            error("The first two input arguments must have the same size.")
    else:
        if _not(isreal(data)):
            theta = (varargin[I[1]] @ M[pi]) / 180
            newdata = data * (cos(theta) + (M[sqrt(-1)] @ sin(theta)))
        else:
            theta = (varargin[I[1]] @ M[pi]) / 180
            m, n = size(data)
            if mod(n, 2) > 0:
                error("The first input arg. must have an EVEN number of columns.")
            eindx = find(mod(M[1:n], 2) == 1)
            nindx = find(mod(M[1:n], 2) == 0)
            east = data[I[:, eindx]]
            north = data[I[:, nindx]]
            newx = (east * cos(theta)) + (north * sin(theta))
            newy = ((-east) * sin(theta)) + (north * cos(theta))
    if length(varargin) > 1:
        varargout[I[1]] = C[newx]
        varargout[I[2]] = C[newy]
    else:
        if _not(isreal(data)):
            varargout[I[1]] = C[newdata]
        else:
            newdata[I[:, eindx]] = copy(newx)
            newdata[I[:, nindx]] = copy(newy)
            varargout[I[1]] = C[newdata]
    return varargout
    return varargout


def cmgpolar(u, v):
    u1 = copy(u)
    v1 = copy(v)
    u1[I[find(isnan(u1))]] = M[[]]
    v1[I[find(isnan(v1))]] = M[[]]
    if all(u1 >= 0) & all(v1 >= 0):
        u = copy(u)
    else:
        u, v = cmguv2spd(u, v)
    v = (M[v] @ pi) / 180
    umax = max(u)
    umin = min(u)
    udiff = umax - umin
    if udiff < 1000:
        umax2 = 100 * ceil(umax / 100)
    if udiff < 100:
        umax2 = 10 * ceil(umax / 10)
    if udiff < 10:
        umax2 = ceil(umax)
    polar(
        M[
            v,
            0,
        ],
        M[
            u,
            umax2,
        ],
        "+",
    )
    view(90, -90)
    return


def cmgpca(u, v, onoff):
    if nargin < 1:
        help(mfilename)
        return
    elif nargin < 2:
        fprintf("\nAt least two input arguments are needed.\n")
        return
    elif nargin < 3:
        onoff = 0
    if any(size(u) - size(v)):
        fprintf("\nSizes of the first two arguments must agree.\n")
    elif size(u, 2) > 1:
        fprintf("\nThe first two arguments cannot be matrices.\n")
    timeplt_figure(12, 1, "anything")
    indx = find(isnan(u) == 0)
    u = u(indx)
    v = v(indx)
    indx = find(isnan(v) == 0)
    u = u(indx)
    v = v(indx)
    mu = mean(u)
    mv = mean(v)
    covar = cov(u, v)
    eigvec, eigval = eig(covar)
    spd, direc = cmguv2spd(mu, mv)
    if eigval(2, 2) > eigval(1, 1):
        cols = M[[2, 1]]
    else:
        cols = M[[1, 2]]
    majr = sqrt(eigval(cols(1), cols(1)))
    junk, majaz = cmguv2spd(eigvec(1, cols(1)), eigvec(2, cols(1)))
    minr = sqrt(eigval(cols(2), cols(2)))
    junk, minaz = cmguv2spd(eigvec(1, cols(2)), eigvec(2, cols(2)))
    set(gcf, "units", "norm")
    hh2[I[1]] = subplot(212)
    set(hh2(1), "position", M[[0.1, 0.05, 0.8, 0.19]])
    hh2[I[2]] = text(
        0.2, 0.8, M[["Speed= ", num2str(spd), ";  Direction= ", num2str(direc)]]
    )
    set(hh2(2), "color", "b")
    hh2[I[3]] = text(
        0.2,
        0.5,
        M[
            [
                "Major axis: Magnitude= ",
                num2str(majr * 2),
                ";  Azimuth= ",
                num2str(majaz),
            ]
        ],
    )
    set(hh2(3), "color", "r")
    hh2[I[4]] = text(
        0.2,
        0.2,
        M[
            [
                "Minor axis: Magnitude= ",
                num2str(minr * 2),
                ";  Azimuth= ",
                num2str(minaz),
            ]
        ],
    )
    set(hh2(4), "color", "r")
    axis("off")
    hand = subplot(211)
    if onoff:
        ticks = ceil(length(u) / 5000)
        plot(u[I[1:ticks:end]], v[I[1:ticks:end]], ".g", "markersize", 1)
    hold("on")
    set(hand, "position", M[[0.1, 0.35, 0.8, 0.6]])
    x1, y1 = cmgspd2uv(majr, majaz)
    x2, y2 = cmgspd2uv(minr, minaz)
    hh1 = plot(
        M[[x1, -x1]],
        M[[y1, -y1]],
        "-r",
        M[[x2, -x2]],
        M[[y2, -y2]],
        "-r",
        "linewidth",
        2,
    )
    theta = ((2 * pi) @ (M[1:64])) / 64
    theta = M[[0, theta]]
    xx = M[majr] @ cos(theta)
    yy = M[minr] @ sin(theta)
    angle = (((-majaz) @ M[pi]) / 180) + (pi / 2)
    xxx = (M[xx] @ cos(angle)) - (M[yy] @ sin(angle))
    yyy = (M[xx] @ sin(angle)) + (M[yy] @ cos(angle))
    hh1[I[3]] = plot(xxx, yyy, "r")
    axis("square")
    if onoff > 0:
        s = min(M[[max(abs(u)), max(abs(v))]])
        s = ceil((0.75 * s) / 10) * 10
    else:
        s = ceil(max(M[[max(M[[x1, y1, x2, y2]]), mu, mv]]) / 10) * 10
    axis(M[[-s, s, -s, s]])
    hh1[I[4]] = plot(
        M[
            0,
            mu,
        ],
        M[
            0,
            mv,
        ],
        "-b",
        "linewidth",
        1,
    )
    xlabel(strrep(inputname(1), "_", "\_"))
    ylabel(strrep(inputname(2), "_", "\_"))
    box("on")
    hold("off")
    return