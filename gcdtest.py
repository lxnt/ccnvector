#!/usr/bin/env python3

import ccnvector as nv

from math import pi as M_PI

WGS84_A = 6378137.0;
WGS84_F = 1/298.257223563;

nv.ellipsoid(WGS84_A, WGS84_F, 400)

d2l = M_PI/180.0

lla = (45.0*d2l, 10.0*d2l)

latb = 45.0*d2l
lonb = 340.0*d2l

na = nv.latlon2n(lla)
if False:
    for lstep in range(45):
        llb = (latb, lonb + lstep*d2l)
        nb = nv.latlon2n(llb)
        
        print("({:.1f} {:.1f})-> ({:.1f} {:.1f}) gcd={:.2f} gcde={:.2f} dcd={:.2f} dcde={:.2f}".format(
            lla[0]/d2l, lla[1]/d2l, llb[0]/d2l, llb[1]/d2l, 
            nv.gcdist(na, nb)/1000, 
            nv.gcdist_ellipsoid(na, nb)/1000,
            nv.cartdist(na, 0.0, nb, 0.0)/1000,
            nv.cartdist_ellipsoid(na, 0.0, nb, 0.0)/1000))
        

print("")
lla = (-45*d2l,-45*d2l)
llb = (45*d2l, 45*d2l)
na=nv.latlon2n(lla)
nb=nv.latlon2n(llb)
print("({:.1f} {:.1f})-> ({:.1f} {:.1f}) gcd={:.2f} gcde={:.2f} dcd={:.2f} dcde={:.2f}".format(
    lla[0]/d2l, lla[1]/d2l, llb[0]/d2l, llb[1]/d2l, 
    nv.gcdist(na, nb)/1000,
    nv.gcdist_ellipsoid(na, nb)/1000,
    nv.cartdist(na, 0.0, nb, 0.0)/1000,
    nv.cartdist_ellipsoid(na, 0.0, nb, 0.0)/1000))

print(nv.n2latlon(na)[0]/d2l, nv.n2latlon(na)[1]/d2l)
print(nv.n2latlon(nb)[0]/d2l, nv.n2latlon(nb)[1]/d2l)
