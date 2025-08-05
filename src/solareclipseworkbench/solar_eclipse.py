# Greg Miller gmiller@gregmiller.net 2023
# https://www.celestialprogramming.com/
# Released as public domain
#
# Functions for computing solar eclipses.
# Algorithms from various sources including:
# - A Manual of Spherical and Practical Astronomy - Chauvenet 1863
# - The Mathematical Theory of Eclipses According to Chauvenet's Transformation of Bessel's Method - Buchanan 1904
# - Explanatory Supplement to the Astronomical Ephemeris 1961
# - Explanatory Supplement to the Astronomical Almanac 1992
# - Elements of Solar Eclipses - Meeus 1989
# - Prediction and Analysis of Solar Eclipse Circumstances - Williams 1971

import math
import csv
import os
from datetime import datetime

rad = math.pi / 180
maxiterations = 20

def solve_quadrant(sin, cos):
    if sin >= 0 and cos >= 0:
        return math.asin(sin)
    if sin < 0 and cos >= 0:
        return math.asin(sin)
    if sin < 0 and cos < 0:
        return -math.acos(cos)
    if sin >= 0 and cos < 0:
        return math.acos(cos)


def get_element_coeffs(date=None):
    """
    Returns eclipse elements for the eclipse closest to the given date.
    If no date is provided, defaults to 2026 eclipse.
    date: str or datetime (YYYY-MM-DD)
    """
    if date is None:
        return get_elements_2026()
    if isinstance(date, str):
        date_obj = datetime.strptime(date, '%Y-%m-%d')
    else:
        date_obj = date

    # Example logic: select by year
    year = date_obj.year
    if year == 2017:
        return get_elements_2017()
    elif year == 2024:
        return get_elements_2024()
    elif year == 1994:
        return get_elements_1994()
    elif year == 1996:
        return get_elements_1996()
    elif year == 2026:
        return get_elements_2026()
    # Default fallback
    return get_elements_2026()


def get_elements_2017():
    elements = {}
    elements['jd'] = 2457987.5
    elements['Δt'] = 69.1
    elements['T0'] = 18

    elements['X0'] = -0.1295710
    elements['X1'] = 0.5406426
    elements['X2'] = -0.0000294
    elements['X3'] = -0.0000081

    elements['Y0'] = 0.4854160
    elements['Y1'] = -0.1416400
    elements['Y2'] = -0.0000905
    elements['Y3'] = 0.0000020

    elements['d0'] = 11.8669596
    elements['d1'] = -0.0136220
    elements['d2'] = -0.0000020
    elements['d3'] = 0.0

    elements['L10'] = 0.5420930
    elements['L11'] = 0.0001241
    elements['L12'] = -0.0000118
    elements['L13'] = 0.0

    elements['L20'] = -0.0040250
    elements['L21'] = 0.0001234
    elements['L22'] = -0.0000117
    elements['L23'] = 0.0

    elements['M0'] = 89.245430
    elements['M1'] = 15.003940
    elements['M2'] = 0.0
    elements['M3'] = 0.0

    elements['tanf1'] = 0.0046222
    elements['tanf2'] = 0.0045992

    return elements


def get_elements_2024_2():
    # year,"month","day","td_ge","dt","luna_num","saros","eclipse_type","gamma","magnitude","lat_ge","lng_ge","lat_dd_ge","lng_dd_ge","sun_alt","sun_azm","path_width","central_duration","duration_secs","cat_no","canon_plate",
    # 2024, 4, 8, "18:18:29", 74.00000000, 300, 139, "T", .34314000, 1.05655000, "25.3N", "104.1W", 25.28945000, -104.12775000, 69.80000000, 149.40000000, 197.50000000, "04m28s", 268.10000000, 9561.00000000, 479.00000000,
    # "julian_date",    "t0","x0","x1","x2","x3","y0","y1","y2","y3","d0","d1","d2","mu0","mu1","mu2","l10","l11","l12","l20","l21","l22","tan_f1","tan_f2","tmin","tmax","etype","PNS","UNS","NCN","nSer","nSeq","nJLE"
    # 2460409.26300000, 18.00000000, -.31824400, .51171160, .00003260, -.00000842, .21976400, .27095890, -.00005950, -.00000466, 7.58620000, .01484400, -.00000200, 89.59122000, 15.00408000, .00000000, .53581400, .00006180, -.00001280, -.01027200, .00006150, -.00001270, .00466830, .00464500, -3.00000000, 3.00000000, 1, 0, 0, 0, 71, 30, 40

    elements = {}
    elements['jd'] = 2460408.5
    elements['Δt'] = 69.1
    elements['T0'] = 18

    elements['X0'] = -0.3182588
    elements['X1'] = 0.5117224
    elements['X2'] = 0.0000330
    elements['X3'] = -0.0000085

    elements['Y0'] = 0.2197690
    elements['Y1'] = 0.2709652
    elements['Y2'] = -0.0000592
    elements['Y3'] = -0.0000047

    elements['d0'] = 7.5861809
    elements['d1'] = 0.0148443
    elements['d2'] = -0.0000017
    elements['d3'] = 0.0

    elements['L10'] = 0.5357259
    elements['L11'] = 0.0000620
    elements['L12'] = -0.0000128
    elements['L13'] = 0.0

    elements['L20'] = -0.0106071
    elements['L21'] = 0.0000617
    elements['L22'] = -0.0000127
    elements['L23'] = 0.0

    elements['M0'] = 89.5912142
    elements['M1'] = 15.0040824
    elements['M2'] = -0.0000008
    elements['M3'] = 0.0

    elements['tanf1'] = 0.0046663
    elements['tanf2'] = 0.0046430

    return elements


def get_elements_2024():
    elements = {}
    elements['jd'] = 2460408.5
    elements['Δt'] = 69.1
    elements['T0'] = 18

    elements['X0'] = -0.3182440
    elements['X1'] = 0.5117116
    elements['X2'] = 0.0000326
    elements['X3'] = -0.0000084

    elements['Y0'] = 0.2197640
    elements['Y1'] = 0.2709589
    elements['Y2'] = -0.0000595
    elements['Y3'] = -0.0000047

    elements['d0'] = 7.5862002
    elements['d1'] = 0.0148440
    elements['d2'] = -0.0000020
    elements['d3'] = 0.0

    elements['L10'] = 0.5358140
    elements['L11'] = 0.0000618
    elements['L12'] = -0.0000128
    elements['L13'] = 0.0

    elements['L20'] = -0.0102720
    elements['L21'] = 0.0000615
    elements['L22'] = -0.0000127
    elements['L23'] = 0.0

    elements['M0'] = 89.591217
    elements['M1'] = 15.004080
    elements['M2'] = 0.0
    elements['M3'] = 0.0

    elements['tanf1'] = 0.0046683
    elements['tanf2'] = 0.0046450

    return elements


def get_elements_1996():
    elements = {}
    elements['Δt'] = 63
    elements['T0'] = 14

    elements['X0'] = 0.296103
    elements['X1'] = 0.5060364
    elements['X2'] = 0.0000145
    elements['X3'] = -0.00000644

    elements['Y0'] = 1.083058
    elements['Y1'] = -0.1515218
    elements['Y2'] = -0.0000102
    elements['Y3'] = 0.00000184

    elements['d0'] = -7.63950
    elements['d1'] = -0.015234
    elements['d2'] = 0.000002

    elements['M0'] = 33.40582
    elements['M1'] = 15.003782

    elements['L10'] = 0.559341
    elements['L11'] = -0.0001067
    elements['L12'] = -0.0000107

    elements['L20'] = 0.013151
    elements['L21'] = -0.0001062
    elements['L22'] = -0.0000106

    elements['tanf1'] = 0.0046865
    elements['tanf2'] = 0.0046632

    return elements


def get_elements_1994():
    elements = {}
    elements['Δt'] = 61
    elements['T0'] = 17

    elements['X0'] = -0.173367
    elements['X1'] = 0.4990629
    elements['X2'] = 0.0000296
    elements['X3'] = -0.00000563

    elements['Y0'] = 0.383484
    elements['Y1'] = 0.0869393
    elements['Y2'] = -0.0001183
    elements['Y3'] = -0.00000092

    elements['d0'] = 17.68613
    elements['d1'] = 0.010642
    elements['d2'] = -0.000004

    elements['M0'] = 75.90923
    elements['M1'] = 15.001621

    elements['L10'] = 0.566906
    elements['L11'] = -0.0000318
    elements['L12'] = -0.0000098

    elements['L20'] = 0.020679
    elements['L21'] = -0.0000317
    elements['L22'] = -0.0000097

    elements['tanf1'] = 0.0046308
    elements['tanf2'] = 0.0046077

    return elements


def get_elements_2026():
    # year,"month","day","td_ge","dt","luna_num","saros","eclipse_type","gamma","magnitude","lat_ge","lng_ge","lat_dd_ge","lng_dd_ge","sun_alt","sun_azm","path_width","central_duration","duration_secs","cat_no","canon_plate",
    # 2026, 8, 12, "17:47:06", 75.40000000, 329, 126, "T", .89774000, 1.03863000, "65.2N", "25.2W", 65.22345000, -25.21619000, 25.80000000, 248.40000000, 293.90000000, "02m18s", 138.20000000, 9566.00000000, 479.00000000,
    # "julian_date",    "t0",        "x0",      "x1",      "x2",       "x3",       "y0",      "y1",       "y2",       "y3",
    # 2461265.24100000, 18.00000000, .47551400, .51892490, -.00007730, -.00000804, .77118300, -.23016800, -.00012460, .00000377
    # "d0"       ,"d1"       ,"d2"       ,"mu0"       ,"mu1"       ,"mu2"     ,"l10"     ,"l11"     ,"l12"      ,"l20"      ,"l21"     ,"l22","tan_f1","tan_f2","tmin","tmax","etype","PNS","UNS","NCN","nSer","nSeq","nJLE"
    # 14.79667000, -.01206500, -.00000300, 88.74779000, 15.00309000, .00000000, .53795500, .00009390, -.00001210, -.00814200, .00009350, -.00001210, .00461410, .00459110, -3.00000000, 3.00000000, 1, 1, 0, 0, 72, 48, 2

    elements = {}
    elements['jd'] = 2461265.241
    elements['Δt'] = 69.184
    elements['T0'] = 18

    elements['X0'] = 0.47551400
    elements['X1'] = 0.51892490
    elements['X2'] = -0.00007730
    elements['X3'] = -0.00000804

    elements['Y0'] = 0.77118300
    elements['Y1'] = -0.23016800
    elements['Y2'] = -0.00012460
    elements['Y3'] = 0.00000377

    elements['d0'] = 14.79667000
    elements['d1'] = -0.01206500
    elements['d2'] = -0.00000300

    elements['L10'] = 0.53795500
    elements['L11'] = 0.00009390
    elements['L12'] = -0.00001210

    # "d0"       ,"d1"       ,"d2"       ,"mu0"       ,"mu1"       ,"mu2"     ,"l10"     ,"l11"     ,"l12"      ,"l20"      ,"l21"     ,"l22"      ,"tan_f1","tan_f2","tmin","tmax","etype","PNS","UNS","NCN","nSer","nSeq","nJLE"
    # 14.79667000, -.01206500, -.00000300, 88.74779000, 15.00309000, .00000000, .53795500, .00009390, -.00001210, -.00814200, .00009350, -.00001210, .00461410, .00459110, -3.00000000, 3.00000000, 1, 1, 0, 0, 72, 48, 2
    elements['L20'] = -0.00814200
    elements['L21'] = 0.00009350
    elements['L22'] = -0.00001210

    elements['M0'] = 88.74779000
    elements['M1'] = 15.00309000
    elements['M2'] = 0.0

    elements['tanf1'] = 0.00461410
    elements['tanf2'] = 0.00459110

    return elements


def get_elements(e, t, phi, lam, height):
    # Meeus - Elements of Solar Eclipses
    o = {}
    o['X'] = e['X0'] + e['X1'] * t + e['X2'] * t * t + e['X3'] * t * t * t
    o['Y'] = e['Y0'] + e['Y1'] * t + e['Y2'] * t * t + e['Y3'] * t * t * t
    o['d'] = e['d0'] + e['d1'] * t + e['d2'] * t * t
    o['M'] = e['M0'] + e['M1'] * t
    o['Xp'] = e['X1'] + 2 * e['X2'] * t + 3 * e['X3'] * t * t
    o['Yp'] = e['Y1'] + 2 * e['Y2'] * t + 3 * e['Y3'] * t * t
    o['Mp'] = e['M1']
    o['L1'] = e['L10'] + e['L11'] * t + e['L12'] * t * t
    o['L2'] = e['L20'] + e['L21'] * t + e['L22'] * t * t
    o['tanf1'] = e['tanf1']
    o['tanf2'] = e['tanf2']

    o['H'] = o['M'] - lam - 0.00417807 * e['Δt']

    o['u1'] = math.atan(0.99664719 * math.tan(phi * rad)) / rad
    o['rho_sin_phip'] = 0.99664719 * math.sin(o['u1'] * rad) + height / 6378140 * math.sin(phi * rad)
    o['rho_cos_phip'] = math.cos(o['u1'] * rad) + height / 6378140 * math.cos(phi * rad)

    o['xi'] = o['rho_cos_phip'] * math.sin(o['H'] * rad)
    o['eta'] = o['rho_sin_phip'] * math.cos(o['d'] * rad) - o['rho_cos_phip'] * math.cos(o['H'] * rad) * math.sin(o['d'] * rad)
    o['zeta'] = o['rho_sin_phip'] * math.sin(o['d'] * rad) + o['rho_cos_phip'] * math.cos(o['H'] * rad) * math.cos(o['d'] * rad)
    o['xi_p'] = 0.01745329 * e['M1'] * o['rho_cos_phip'] * math.cos(o['H'] * rad)
    o['eta_p'] = 0.01745329 * (e['M1'] * o['xi'] * math.sin(o['d'] * rad) - o['zeta'] * e['d1'])
    o['L1p'] = o['L1'] - o['zeta'] * e['tanf1']
    o['L2p'] = o['L2'] - o['zeta'] * e['tanf2']

    o['u'] = o['X'] - o['xi']
    o['v'] = o['Y'] - o['eta']
    o['a'] = o['Xp'] - o['xi_p']
    o['b'] = o['Yp'] - o['eta_p']
    o['n'] = math.sqrt(o['a'] * o['a'] + o['b'] * o['b'])

    return o


def get_local_circumstances(phi, lam, height, date = None):
    # Meeus - Elements of Solar Eclipses
    e = get_element_coeffs(date)
    lam = -lam

    t = 0
    tau_m = 10000
    iterations = 0
    o = None
    while abs(tau_m) > 0.00001 and iterations < maxiterations:
        o = get_elements(e, t, phi, lam, height)
        tau_m = - (o['u'] * o['a'] + o['v'] * o['b']) / (o['n'] * o['n'])
        t = t + tau_m
        iterations += 1

    m = math.sqrt(o['u'] * o['u'] + o['v'] * o['v'])
    G = (o['L1p'] - m) / (o['L1p'] + o['L2p'])
    A = (o['L1p'] - o['L2p']) / (o['L1p'] + o['L2p'])

    Pm = math.atan2(o['u'] / o['v'], 1) / rad

    sinh = math.sin(o['d'] * rad) * math.sin(phi * rad) + math.cos(o['d'] * rad) * math.cos(phi * rad) * math.cos(o['H'] * rad)
    h = math.asin(sinh) / rad
    q = math.asin(math.cos(phi * rad) * math.sin(o['H'] * rad) / math.cos(h * rad))

    S = (o['a'] * o['v'] - o['u'] * o['b']) / (o['n'] * o['L1p'])
    tau = o['L1p'] / o['n'] * math.sqrt(max(0, 1 - S * S))

    first_contact = t - tau
    last_contact = t + tau

    for _ in range(10):
        fco = get_elements(e, first_contact, phi, lam, height)
        S = (fco['a'] * fco['v'] - fco['u'] * fco['b']) / (fco['n'] * fco['L1p'])
        tau_f = - (fco['u'] * fco['a'] + fco['v'] * fco['b']) / (fco['n'] * fco['n']) - fco['L1p'] / fco['n'] * math.sqrt(max(0, 1 - S * S))
        first_contact = first_contact + tau_f

    for _ in range(10):
        fco = get_elements(e, last_contact, phi, lam, height)
        S = (fco['a'] * fco['v'] - fco['u'] * fco['b']) / (fco['n'] * fco['L1p'])
        tau_f = - (fco['u'] * fco['a'] + fco['v'] * fco['b']) / (fco['n'] * fco['n']) + fco['L1p'] / fco['n'] * math.sqrt(max(0, 1 - S * S))
        last_contact = last_contact + tau_f

    # Interior contacts
    S = (o['a'] * o['v'] - o['u'] * o['b']) / (o['n'] * o['L2p'])
    tau = o['L2p'] / o['n'] * math.sqrt(max(0, 1 - S * S))

    third_contact = t - tau
    second_contact = t + tau

    for _ in range(10):
        fco = get_elements(e, third_contact, phi, lam, height)
        S = (fco['a'] * fco['v'] - fco['u'] * fco['b']) / (fco['n'] * fco['L2p'])
        tau_f = - (fco['u'] * fco['a'] + fco['v'] * fco['b']) / (fco['n'] * fco['n']) - fco['L2p'] / fco['n'] * math.sqrt(max(0, 1 - S * S))
        third_contact = third_contact + tau_f

    for _ in range(10):
        fco = get_elements(e, second_contact, phi, lam, height)
        S = (fco['a'] * fco['v'] - fco['u'] * fco['b']) / (fco['n'] * fco['L2p'])
        tau_f = - (fco['u'] * fco['a'] + fco['v'] * fco['b']) / (fco['n'] * fco['n']) + fco['L2p'] / fco['n'] * math.sqrt(max(0, 1 - S * S))
        second_contact = second_contact + tau_f

    UTFirstContact = e['T0'] + first_contact - e['Δt'] / 60 / 60
    UTSecondContact = e['T0'] + second_contact - e['Δt'] / 60 / 60
    UTThirdContact = e['T0'] + third_contact - e['Δt'] / 60 / 60
    UTLastContact = e['T0'] + last_contact - e['Δt'] / 60 / 60
    UTMaximum = e['T0'] + t - e['Δt'] / 60 / 60

    UT = e['T0'] + t

    return {
        'jd': e['jd'], 't': t, 'UT': UT, 'mag': G,
        'UTMaximum': UTMaximum,
        'UTFirstContact': UTFirstContact,
        'UTSecondContact': UTSecondContact,
        'UTThirdContact': UTThirdContact,
        'UTLastContact': UTLastContact,
        'h': h, 'm': m, 'elements': o
    }


def compute_central_lat_lon_for_time(e, UTC):
    t = UTC - e['T0']
    X = e['X0'] + e['X1'] * t + e['X2'] * t * t + e['X3'] * t * t * t
    Y = e['Y0'] + e['Y1'] * t + e['Y2'] * t * t + e['Y3'] * t * t * t
    d = e['d0'] + e['d1'] * t + e['d2'] * t * t
    M = e['M0'] + e['M1'] * t
    L1 = e['L10'] + e['L11'] * t + e['L12'] * t * t
    L2 = e['L20'] + e['L21'] * t + e['L22'] * t * t
    Xp = e['X1'] + 2 * e['X2'] * t + 3 * e['X3'] * t * t
    Yp = e['Y1'] + 2 * e['Y2'] * t + 3 * e['Y3'] * t * t
    dtemp = d * rad
    omega = 1 / math.sqrt(1 - 0.006694385 * math.cos(d * rad) * math.cos(d * rad))
    p = e['M1'] / 57.2957795
    b = Yp - p * X * math.sin(d * rad)
    c = Xp + p * Y * math.sin(d * rad)
    y1 = omega * Y
    b1 = omega * math.sin(d * rad)
    b2 = 0.99664719 * omega * math.cos(d * rad)
    B = math.sqrt(max(0, 1 - X * X - y1 * y1))
    Phi1 = math.asin(B * b1 + y1 * b2)
    sinH = X / math.cos(Phi1)
    cosH = (B * b2 - y1 * b1) / math.cos(Phi1)
    H = solve_quadrant(sinH, cosH) / rad
    Phi = math.atan(1.00336409 * math.tan(Phi1)) / rad
    lam = M - H - 0.00417807 * e['Δt']
    L1p = L1 - B * e['tanf1']
    L2p = L2 - B * e['tanf2']
    a = c - p * B * math.cos(d * rad)
    n = math.sqrt(a * a + b * b)
    duration = 7200 * L2p / n
    sinh = math.sin(d * rad) * math.sin(Phi * rad) + math.cos(d * rad) * math.cos(Phi * rad) * math.cos(H * rad)
    h = math.asin(sinh) / rad
    K2 = B * B + ((X * a + Y * b) * (X * a + Y * b)) / (n * n)
    K = math.sqrt(K2)
    width = 12756 * abs(L2p) / K
    A = (L1p - L2p) / (L1p + L2p)
    return {'lat': Phi, 'lon': -lam, 'magnitude': A, 'duration': duration, 'width': width}


def compute_extremes(e, t0):
    t = t0
    r = {}
    r['t'] = t
    r['X'] = e['X0'] + e['X1'] * t + e['X2'] * t * t + e['X3'] * t * t * t
    r['Y'] = e['Y0'] + e['Y1'] * t + e['Y2'] * t * t + e['Y3'] * t * t * t
    r['d'] = e['d0'] + e['d1'] * t + e['d2'] * t * t
    r['M'] = e['M0'] + e['M1'] * t
    r['L1'] = e['L10'] + e['L11'] * t + e['L12'] * t * t
    r['L2'] = e['L20'] + e['L21'] * t + e['L22'] * t * t
    r['Xp'] = e['X1'] + 2 * e['X2'] * t + 3 * e['X3'] * t * t
    r['Yp'] = e['Y1'] + 2 * e['Y2'] * t + 3 * e['Y3'] * t * t
    r['omega'] = 1 / math.sqrt(1 - 0.006694385 * math.cos(r['d'] * rad) * math.cos(r['d'] * rad))
    r['p'] = e['M1'] / 57.2957795
    r['b'] = r['Yp'] - r['p'] * r['X'] * math.sin(r['d'] * rad)
    r['c'] = r['Xp'] + r['p'] * r['Y'] * math.sin(r['d'] * rad)
    r['y1'] = r['omega'] * r['Y']
    r['b1'] = r['omega'] * math.sin(r['d'] * rad)
    r['b2'] = 0.99664719 * r['omega'] * math.cos(r['d'] * rad)
    temp = 1 - r['X'] * r['X'] - r['y1'] * r['y1']
    if temp < 0:
        temp = 0
    r['B'] = math.sqrt(temp)
    r['Phi1'] = math.asin(r['B'] * r['b1'] + r['y1'] * r['b2'])
    sinH = r['X'] / math.cos(r['Phi1'])
    cosH = (r['B'] * r['b2'] - r['y1'] * r['b1']) / math.cos(r['Phi1'])
    r['H'] = solve_quadrant(sinH, cosH) / rad
    r['Phi'] = math.atan(1.00336409 * math.tan(r['Phi1'])) / rad
    r['lam'] = -(r['M'] - r['H'] - 0.00417807 * e['Δt'])
    r['L1p'] = r['L1'] - r['B'] * e['tanf1']
    r['L2p'] = r['L2'] - r['B'] * e['tanf2']
    r['a'] = r['c'] - r['p'] * r['B'] * math.cos(r['d'] * rad)
    r['n'] = math.sqrt(r['a'] * r['a'] + r['b'] * r['b'])
    r['duration'] = 7200 * r['L2p'] / r['n']
    r['sinh'] = math.sin(r['d'] * rad) * math.sin(r['Phi'] * rad) + math.cos(r['d'] * rad) * math.cos(r['Phi'] * rad) * math.cos(r['H'] * rad)
    r['h'] = math.asin(r['sinh']) / rad
    r['K2'] = r['B'] * r['B'] + ((r['X'] * r['a'] + r['Y'] * r['b']) * (r['X'] * r['a'] + r['Y'] * r['b'])) / (r['n'] * r['n'])
    r['K'] = math.sqrt(r['K2'])
    r['width'] = 12756 * abs(r['L2p']) / r['K']
    r['A'] = (r['L1p'] - r['L2p']) / (r['L1p'] + r['L2p'])
    return r


def compute_estimate(e):
    omega = 1 / math.sqrt(1 - 0.006694385 * math.cos(e['d0'] * rad) * math.cos(e['d0'] * rad))
    u = e['X0']
    a = e['X1']
    v = omega * e['Y0']
    b = omega * e['Y1']
    n = math.sqrt(a * a + b * b)
    S = (a * v - u * b) / n
    tau = -(u * a + v * b) / (n * n)
    tau1 = tau - math.sqrt(1 - S * S) / n
    tau2 = tau + math.sqrt(1 - S * S) / n
    return {'tau1': tau1, 'tau2': tau2}


def refine_estimate(e, t):
    X = e['X0'] + e['X1'] * t + e['X2'] * t * t + e['X3'] * t * t * t
    Y = e['Y0'] + e['Y1'] * t + e['Y2'] * t * t + e['Y3'] * t * t * t
    d = e['d0'] + e['d1'] * t + e['d2'] * t * t
    Xp = e['X1'] + 2 * e['X2'] * t + 3 * e['X3'] * t * t
    Yp = e['Y1'] + 2 * e['Y2'] * t + 3 * e['Y3'] * t * t
    omega = 1 / math.sqrt(1 - 0.006694385 * math.cos(e['d0'] * rad) * math.cos(e['d0'] * rad))
    u = X
    a = Xp
    v = omega * Y
    b = omega * Yp
    n = math.sqrt(a * a + b * b)
    S = (a * v - u * b) / n
    tau = -(u * a + v * b) / (n * n)
    tau1 = tau - math.sqrt(1 - S * S) / n
    tau2 = tau + math.sqrt(1 - S * S) / n
    return {'tau1': tau1, 'tau2': tau2}


def get_extreme_points(e):
    est = compute_estimate(e)
    est2 = refine_estimate(e, est['tau1'])
    est3 = refine_estimate(e, est['tau2'])
    begin_circumstances = compute_extremes(e, est['tau1'] + est2['tau1'])
    end_circumstances = compute_extremes(e, est['tau2'] + est3['tau2'])
    return {'begin': begin_circumstances, 'end': end_circumstances}


def compute_rise_set_point(be, gamma):
    eta = math.cos(gamma)
    xi = math.sin(gamma)
    sind = math.sin(be['d'] * rad)
    cosd = math.cos(be['d'] * rad)
    cosphi1sintheta = xi
    cosphi1costheta = -eta * sind
    theta = math.atan2(cosphi1sintheta, cosphi1costheta)
    lam = -(be['H'] * rad - theta)
    sinphi = eta * cosd
    phi = math.asin(sinphi)
    return {'lat': phi / rad, 'lon': lam / rad if lam <= math.pi else lam / rad - 360}

def compute_rise_set_points(be):
    m = math.sqrt(be['X'] * be['X'] + be['Y'] * be['Y'])
    M = math.atan2(be['X'], be['Y'])
    cos_gamma_M = (m * m + 1 - be['L1'] * be['L1']) / (2 * m)
    gamma_M = math.acos(cos_gamma_M)
    gamma_M2 = 2 * math.pi - gamma_M
    gamma1 = gamma_M + M
    gamma2 = gamma_M2 + M
    ll1 = compute_rise_set_point(be, gamma1)
    ll2 = compute_rise_set_point(be, gamma2)
    return [ll1, ll2]

def get_rise_set_curves():
    t = -3
    e = get_element_coeffs()
    list1 = []
    list2 = []
    list3 = []
    list4 = []
    maxlat = -100
    minlat = 100
    nStart = None
    sStart = None
    while t < .5:
        be = get_elements(e, t, 0, 0, 0)
        p = compute_rise_set_points(be)
        if not math.isnan(p[0]['lat']):
            list1.append(p[0])
            if p[0]['lat'] > maxlat:
                maxlat = p[0]['lat']
                nStart = p[0]['lon']
            if p[0]['lat'] < minlat:
                minlat = p[0]['lat']
                sStart = p[0]['lon']
        if not math.isnan(p[1]['lat']):
            list2.insert(0, p[1])
            if p[1]['lat'] < minlat:
                minlat = p[1]['lat']
                sStart = p[1]['lon']
        t += .01
    t = .5
    maxlat = -100
    minlat = 100
    nEnd = None
    sEnd = None
    while t < 4:
        be = get_elements(e, t, 0, 0, 0)
        p = compute_rise_set_points(be)
        if not math.isnan(p[0]['lat']):
            list3.append(p[0])
            if p[0]['lat'] > maxlat:
                maxlat = p[0]['lat']
                nEnd = p[0]['lon']
            if p[0]['lat'] < minlat:
                minlat = p[0]['lat']
                sEnd = p[0]['lon']
        if not math.isnan(p[1]['lat']):
            list4.insert(0, p[1])
            if p[1]['lat'] > maxlat:
                maxlat = p[1]['lat']
                nEnd = p[1]['lon']
            if p[1]['lat'] < minlat:
                minlat = p[1]['lat']
                sEnd = p[1]['lon']
        t += .01
    return {'setting': list1 + list2, 'rising': list3 + list4, 'nStart': nStart, 'sStart': sStart, 'nEnd': nEnd, 'sEnd': sEnd}

def get_limits_by_longitude_as_list(e, northsouth, G, start_lon, end_lon):
    eq_points = []
    polar_points = []
    step = 0.01
    i = start_lon
    while i <= end_lon:
        eq = get_limits_for_longitude(e, i, northsouth, G, 0)
        polar = get_limits_for_longitude(e, i, northsouth, G, 89.9 * (1 if e['Y0'] >= 0 else -1))
        if polar is not None and eq is not None:
            if abs(eq['lat'] - polar['lat']) < 0.1:
                eq_points.append(eq)
            else:
                eq_points.append(eq)
                polar_points.append(polar)
        else:
            eq_points.append(None)
            polar_points.append(None)
        i += step
    return [eq_points, polar_points]

def get_limits_for_longitude(e, lam, northsouth, G, start_phi):
    t = 0
    phi = start_phi
    i = 0
    delta_phi = 1000
    tau = 1000
    while (abs(tau) > 0.0001 or abs(delta_phi) > 0.0001) and i < 20:
        X = e['X0'] + e['X1'] * t + e['X2'] * t * t + e['X3'] * t * t * t
        Y = e['Y0'] + e['Y1'] * t + e['Y2'] * t * t + e['Y3'] * t * t * t
        d = e['d0'] + e['d1'] * t + e['d2'] * t * t
        M = e['M0'] + e['M1'] * t
        Xp = e['X1'] + 2 * e['X2'] * t + 3 * e['X3'] * t * t
        Yp = e['Y1'] + 2 * e['Y2'] * t + 3 * e['Y3'] * t * t
        L1 = e['L10'] + e['L11'] * t + e['L12'] * t * t
        L2 = e['L20'] + e['L21'] * t + e['L22'] * t * t
        H = M + lam - 0.00417807 * e['Δt']
        height = 0
        u1 = math.atan(0.99664719 * math.tan(phi * rad)) / rad
        rho_sin_phip = 0.99664719 * math.sin(u1 * rad) + height / 6378140 * math.sin(phi * rad)
        rho_cos_phip = math.cos(u1 * rad) + height / 6378140 * math.cos(phi * rad)
        xi = rho_cos_phip * math.sin(H * rad)
        eta = rho_sin_phip * math.cos(d * rad) - rho_cos_phip * math.cos(H * rad) * math.sin(d * rad)
        zeta = rho_sin_phip * math.sin(d * rad) + rho_cos_phip * math.cos(H * rad) * math.cos(d * rad)
        xi_p = 0.01745329 * e['M1'] * rho_cos_phip * math.cos(H * rad)
        eta_p = 0.01745329 * (e['M1'] * xi * math.sin(d * rad) - zeta * e['d1'])
        L1p = L1 - zeta * e['tanf1']
        L2p = L2 - zeta * e['tanf2']
        u = X - xi
        v = Y - eta
        a = Xp - xi_p
        b = Yp - eta_p
        n = math.sqrt(a * a + b * b)
        tau = - (u * a + v * b) / (n * n)
        W = (v * a - u * b) / n
        Q = ((b * math.sin(H * rad) * rho_sin_phip + a * (math.cos(H * rad) * math.sin(d * rad) * rho_sin_phip + math.cos(d * rad) * rho_cos_phip)))/(57.29578 * n)
        E = L1p - G * (L1p + L2p)
        delta_phi = (W + northsouth * abs(E)) / Q
        t = t + tau
        phi = phi + delta_phi
        i += 1
    if abs(tau) > 0.0001 or abs(delta_phi) > 0.0001:
        return None
    UT = e['T0'] + t
    phi = (90 + phi) % 180
    if phi < 0:
        phi += 180
    phi -= 90
    return {'t': t, 'lat': phi, 'lon': lam}

def compute_outline_point(be, Q, umbra):
    # The Explanatory Supplement to the Astronomical Ephemeris 1961
    e = math.sqrt(0.00672267)
    sind = math.sin(be['d'] * rad)
    cosd = math.cos(be['d'] * rad)
    rho1 = math.sqrt(1 - e * e * cosd * cosd)
    rho2 = math.sqrt(1 - e * e * sind * sind)
    sind1 = sind / rho1
    cosd1 = math.sqrt(1 - e * e) * cosd / rho1
    sind1d2 = e * e * sind * cosd / (rho1 * rho2)
    cosd1d2 = math.sqrt(1 - e * e) / (rho1 * rho2)
    Q *= rad
    sinQ = math.sin(Q)
    cosQ = math.cos(Q)
    tanf = be['tanf1']
    l = be['L1']
    if umbra:
        l = be['L2']
        tanf = be['tanf2']
    xi = be['X'] - l * sinQ
    eta = (be['Y'] - l * cosQ) / rho1
    zeta1 = math.sqrt(1 - xi * xi - eta * eta)
    zeta = rho2 * (zeta1 * cosd1d2 - eta * sind1d2)
    L = l - zeta * tanf
    xi = be['X'] - L * sinQ
    eta = (be['Y'] - L * cosQ) / rho1
    zeta1 = math.sqrt(1 - xi * xi - eta * eta)
    zeta = rho2 * (zeta1 * cosd1d2 - eta * sind1d2)
    L = l - zeta * tanf
    xi = be['X'] - L * sinQ
    eta = (be['Y'] - L * cosQ) / rho1
    zeta1 = math.sqrt(1 - xi * xi - eta * eta)
    cosphi1sintheta = xi
    cosphi1costheta = zeta1 * cosd1 - eta * sind1
    theta = math.atan2(cosphi1sintheta, cosphi1costheta) / rad
    lam = be['H'] - theta
    phi1 = math.asin(eta * cosd1 + zeta1 * sind1)
    phi = math.atan((1 / math.sqrt(1 - e * e)) * math.tan(phi1)) / rad
    if phi > 90:
        phi -= 180
    if phi < -90:
        phi += 180
    if lam > 180:
        lam -= 360
    return {'lat': phi, 'lon': -lam}

def propper_angle(d):
    t = d
    if t < 0:
        t += 360
    if t >= 360:
        t -= 360
    return t

def get_solar_eclipses(number_of_eclipses=None, start_date=None):
    """
    Reads eclipse_besselian.csv and returns a list of solar eclipses.
    Optional parameters:
        number_of_eclipses: int, limits the number of returned eclipses
        start_date: str or datetime, filters eclipses after this date (YYYY-MM-DD)
    """
    eclipses = []
    csv_path = os.path.join(os.path.dirname(__file__), 'eclipse_besselian.csv')
    with open(csv_path, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            eclipse_date = datetime(int(row['year']), int(row['month']), int(row['day']))
            if start_date:
                if isinstance(start_date, str):
                    start_dt = datetime.strptime(start_date, '%Y-%m-%d')
                else:
                    start_dt = start_date
                if eclipse_date < start_dt:
                    continue
            eclipses.append(row)
            if number_of_eclipses and len(eclipses) >= number_of_eclipses:
                break
    return eclipses

def get_outline_curve_q_range(be, l):
    msq = be['X'] * be['X'] + be['Y'] * be['Y']
    m = math.sqrt(msq)
    M = math.atan2(be['X'], be['Y'])
    denom = 2 * l * m
    numer = m * m + l * l - 1
    cosQM = numer / denom
    try:
        Q1 = propper_angle((math.acos(cosQM) + M) / rad)
        Q2 = propper_angle((-math.acos(cosQM) + M) / rad)
    except ValueError:
        Q2 = 0
        Q1 = 360
    if Q1 < Q2:
        Q1 += 360
    return {'start': Q2, 'end': Q1}
