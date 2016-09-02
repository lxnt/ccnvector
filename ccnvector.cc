/*

n_vector_library python interface, minimal.
plus some helpers.

WGS84 is default.

Eigen::Vector3d lat_long2n_E( const double& latitude, const double& longitude );

    convert lan/lon to nvector

    (3-tuple) lat_long2n_E(lat, lon)

LatLongPosition n_E2lat_long( const Eigen::Vector3d& n_E );

    convert nvector to lan/lon

    (2-tuple) n_E2lat_long((3-tuple) n_E)

Eigen::Vector3d n_EA_E_and_n_EB_E2p_AB_E( const Eigen::Vector3d& n_EA_E,
                                          Eigen::Vector3d& n_EB_E,
                                          const double& z_EA = 0,
                                          const double& z_EB = 0,
                                          const double& a = WGS84_A,
                                          const double& f = WGS84_F );

    From two positions A and B, finds the delta position.

    (3-tuple) n_EA_E_and_n_EB_E2p_AB_E ((3-tuple) n_EA_E,
                                        (3-tuple) n_EB_E,
                                        z_EA
                                        z_EB);


Eigen::Vector3d n_EA_E_and_p_AB_E2n_EB_E( const Eigen::Vector3d& n_EA_E,
                                          const Eigen::Vector3d& p_AB_E,
                                          const double& z_EA = 0,
                                          const double& a = WGS84_A,
                                          const double& f = WGS84_F );

    From position A and delta, finds position B. (without Z)

    (3-tuple) n_EA_E_and_p_AB_E2n_EB_E((3-tuple) n_EA_E, (3-tuple) p_AB_E, z_EA);

NVectorPosition n_EA_E_and_p_AB_E2n_EB_E_and_z_EB( const Eigen::Vector3d& n_EA_E,
                                                   const Eigen::Vector3d& p_AB_E,
                                                   const double& z_EA = 0,
                                                   const double& a = WGS84_A,
                                                   const double& f = WGS84_F );

    From position A and delta, finds position B. (with Z)

    (4-tuple) n_EA_E_and_p_AB_E2n_EB_E_and_z_EB((3-tuple) n_EA_E, (3-tuple) p_AB_E, z_EA);


========

Eigen::Vector3d n_EB_E2p_EB_E( const Eigen::Vector3d& n_EB_E,
                               const double& z_EB = 0,
                               const double& a = WGS84_A,
                               const double& f = WGS84_F );

    (3-tuple) n_EB_E2p_EB_E((3-tuple) n_EB_E, z_EB)

        Converts n-vector to Cartesian position vector in meters.

        not needed.

Eigen::Vector3d p_EB_E2n_EB_E( const Eigen::Vector3d& p_EB_E, const double& a = WGS84_A, const double& f = WGS84_F );

    (3-tuple) n_EB_E2p_EB_E((3-tuple) p_EB_E)

        Converts Cartesian position vector in meters to n-vector. (without depth)

        not needed.


NVectorPosition p_EB_E2n_EB_E_and_z_EB( const Eigen::Vector3d& p_EB_E, const double& a = WGS84_A, const double& f = WGS84_F );

    (4-tuple) p_EB_E2n_EB_E_and_z_EB((3-tuple) p_EB_E);

        Converts Cartesian position vector in meters to n-vector. (with depth)

        not needed


and other matrix stuff also not needed (now)

*/

#include "Python.h"
#include "n_vector_library/n_vector_library.h"

using namespace n_vector_library;

// a -eq radius
// c -polar radius
// flattening = (a-c)/a
// a*flattening = a - c
// a*flattening - a = -c
// c= a - a * flattening
// c = a(1 -flattening)
// c = a(1 - 1/etm)
// c= a - a/etm
// Rmean = (2a + c)/3


static double ELL_A = WGS84_A; // metres
static double ELL_F = WGS84_F;
static double ELL_C = ELL_A - ELL_A * ELL_F;
static double ELL_RMEAN = (2.0 * ELL_A + ELL_C)/3.0;
static double ELL_FLATSTEP = 4000.0;  // assuming 4 km is flat enough

static PyObject *
ccnv_ellipsoid(PyObject *self, PyObject *args)
{
    double f = nan("");
    double fstep = nan("");

    if (!PyArg_ParseTuple(args, "d|dd:ellipsoid", &ELL_A, &f, &fstep))
        return NULL;

    if (not isnan(f))
        ELL_F = f;

    if (not isnan(fstep))
        ELL_FLATSTEP = fstep;

    ELL_C = ELL_A - ELL_A * ELL_F;
    ELL_RMEAN = (2.0 * ELL_A + ELL_C)/3.0;
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *
ccnv_lat_long2n_E(PyObject *self, PyObject *args)
{
    double latitude, longitude;

    if (!PyArg_ParseTuple(args, "(dd):latlon2n", &latitude, &longitude))
        return NULL;

    Eigen::Vector3d ev = lat_long2n_E( latitude, longitude );

    return Py_BuildValue("ddd", ev[0], ev[1], ev[2]);
}

static PyObject *
ccnv_n_E2lat_long(PyObject *self, PyObject *args)
{
    double e0, e1, e2;

    if (!PyArg_ParseTuple(args, "(ddd):n2latlon", &e0, &e1, &e2))
        return NULL;
    Eigen::Vector3d ev;
    ev[0] = e0;
    ev[1] = e1;
    ev[2] = e2;
    LatLongPosition llp = n_E2lat_long( ev );

    return Py_BuildValue("dd", llp.latitude, llp.longitude);
}

static PyObject * // AB2D(A, B, zA, zB) ->p_ED
ccnv_n_EA_E_and_n_EB_E2p_AB_E(PyObject *self, PyObject *args)
{
    Eigen::Vector3d ea;
    Eigen::Vector3d eb;
    double a0, a1, a2, b0, b1, b2, az, bz;

    if (!PyArg_ParseTuple(args, "(ddd)(ddd)dd", &a0, &a1, &a2, &b0, &b1, &b2, &az, &bz))
        return NULL;

    ea[0] = a0; ea[1] = a1; ea[2] = a2;
    eb[0] = b0; eb[1] = b1; eb[2] = b2;

    Eigen::Vector3d ev = n_EA_E_and_n_EB_E2p_AB_E( ea, eb, az, bz, ELL_A, ELL_F );

    return Py_BuildValue("ddd", ev[0], ev[1], ev[2]);
}

static PyObject * // AD2B(A, D, zA) -> n_EB
ccnv_n_EA_E_and_p_AB_E2n_EB_E(PyObject *self, PyObject *args)
{
    Eigen::Vector3d ea;
    Eigen::Vector3d eb;
    double a0, a1, a2, b0, b1, b2, az;

    if (!PyArg_ParseTuple(args, "(ddd)(ddd)d", &a0, &a1, &a2, &b0, &b1, &b2, &az))
        return NULL;

    ea[0] = a0; ea[1] = a1; ea[2] = a2;
    eb[0] = b0; eb[1] = b1; eb[2] = b2;

    Eigen::Vector3d ev = n_EA_E_and_p_AB_E2n_EB_E( ea, eb, az, ELL_A, ELL_F );

    return Py_BuildValue("ddd", ev[0], ev[1], ev[2]);
}

static PyObject * // AD2Bz(A, D, zA) -> n_EB, zB
ccnv_n_EA_E_and_p_AB_E2n_EB_E_and_z_EB(PyObject *self, PyObject *args)
{
    Eigen::Vector3d ea;
    Eigen::Vector3d eb;
    double a0, a1, a2, b0, b1, b2, az;

    if (!PyArg_ParseTuple(args, "(ddd)(ddd)d", &a0, &a1, &a2, &b0, &b1, &b2, &az))
        return NULL;

    ea[0] = a0; ea[1] = a1; ea[2] = a2;
    eb[0] = b0; eb[1] = b1; eb[2] = b2;

    NVectorPosition nv = n_EA_E_and_p_AB_E2n_EB_E_and_z_EB( ea, eb, az, ELL_A, ELL_F );

    return Py_BuildValue("(ddd)d", nv.n_EB_E[0], nv.n_EB_E[1], nv.n_EB_E[2], nv.z_EB);
}

static PyObject * // sizevec((x, y, z), l) ->(x, y, z), normalize, then multiply by l
ccnv_sizevec(PyObject *self, PyObject *args)
{
    Eigen::Vector3d ev;
    double a0, a1, a2, l;

    if (!PyArg_ParseTuple(args, "(ddd)d:sizevec", &a0, &a1, &a2, &l))
        return NULL;

    ev[0] = a0; ev[1] = a1; ev[2] = a2;

    ev.normalize();
    ev *= l;

    return Py_BuildValue("ddd", ev[0], ev[1], ev[2]);
}

static PyObject * // advance2d((na0,na1,na2), (nb0,nb1,nb2), l) -> ((nc0, nc1, nc2), (lat, lon)) : advance na in direction of nb for l units over gc line.
ccnv_advance2d(PyObject *self, PyObject *args)
{
    Eigen::Vector3d ea;
    Eigen::Vector3d eb;
    double a0, a1, a2, b0, b1, b2, l;

    if (!PyArg_ParseTuple(args, "(ddd)(ddd)d:advance2d", &a0, &a1, &a2, &b0, &b1, &b2, &l))
        return NULL;

    ea[0] = a0; ea[1] = a1; ea[2] = a2;
    eb[0] = b0; eb[1] = b1; eb[2] = b2;

    Eigen::Vector3d ed = n_EA_E_and_n_EB_E2p_AB_E( ea, eb, 0, 0, ELL_A, ELL_F );

    ed.normalize();
    ed *= l;

    Eigen::Vector3d ena = n_EA_E_and_p_AB_E2n_EB_E( ea, ed, 0, ELL_A, ELL_F );

    LatLongPosition llp = n_E2lat_long( ena );

    return Py_BuildValue("(ddd)(dd)", ena[0], ena[1], ena[2], llp.latitude, llp.longitude);
}

static PyObject * // gcdist((na0,na1,na2), (nb0,nb1,nb2)) -> distance : return great circle distance
ccnv_gcdist(PyObject *self, PyObject *args)
{
    Eigen::Vector3d ea;
    Eigen::Vector3d eb;
    double a0, a1, a2, b0, b1, b2;

    if (!PyArg_ParseTuple(args, "(ddd)(ddd):gcdist", &a0, &a1, &a2, &b0, &b1, &b2))
        return NULL;

    ea[0] = a0; ea[1] = a1; ea[2] = a2;
    eb[0] = b0; eb[1] = b1; eb[2] = b2;

    double rv = atan2(ea.cross(eb).norm(), ea.dot(eb)) * ELL_RMEAN;

    return Py_BuildValue("d", rv);
}
static void print_ll(const Eigen::Vector3d &n)
{

    LatLongPosition llp = n_E2lat_long( n );

    std::cout<<"("<<llp.latitude<<","<<llp.longitude<<")";
}
static PyObject * // gcdist_ellipsoid((na0,na1,na2), (nb0,nb1,nb2)[, r]) -> distance : return great circle
ccnv_gcdist_ellipsoid(PyObject *self, PyObject *args)
{
    Eigen::Vector3d ea;
    Eigen::Vector3d eb;
    double a0, a1, a2, b0, b1, b2;

    if (!PyArg_ParseTuple(args, "(ddd)(ddd):gcdist", &a0, &a1, &a2, &b0, &b1, &b2))
        return NULL;

    /* The Inverse Problem

        Either pull in common solution, like in the geographiclib
        which is iterative and a ton of code, not to mention the dependency.
        That's what they did in py-nvector, lazy cheats.

        Or iterate here using advance2d in small steps
        like quarter degree or something, and summing the distance 'travelled'.
    */
    ea[0] = a0; ea[1] = a1; ea[2] = a2;
    eb[0] = b0; eb[1] = b1; eb[2] = b2;

    // if distance isn't significant, just use sphere.
    double d_sph = atan2(ea.cross(eb).norm(), ea.dot(eb)) * ELL_RMEAN;
    if (d_sph < 2 * ELL_FLATSTEP)
        return Py_BuildValue("d", d_sph);

    std::cout<<" d_sph: "<<(d_sph/1000)<<" ELL_FLATSTEP "<<(ELL_FLATSTEP/1000)<<std::endl;

    double rv = 0.0;
    Eigen::Vector3d ethis = ea;
    Eigen::Vector3d enext = ea;
    Eigen::Vector3d edelta;
    double rsum = 0.0;
    //std::cout<<steps<<" steps of "<<step<<std::endl;
    while(true)
    {
        edelta = n_EA_E_and_n_EB_E2p_AB_E( ethis, eb, 0, 0, ELL_A, ELL_F );
        //std::cout<<edelta.norm()<<std::endl;
        edelta.normalize();
        //std::cout<<edelta.norm()<<std::endl;
        edelta *= ELL_FLATSTEP;
        std::cout<<" ed "<<(edelta.norm()/1000)<<" vs ";
        rsum += edelta.norm();
        ethis = enext;
        enext = n_EA_E_and_p_AB_E2n_EB_E( ethis, edelta, 0, ELL_A, ELL_F );
        edelta = enext - ethis;
        rv += edelta.norm() * ELL_RMEAN;
        std::cout<<(edelta.norm() * ELL_RMEAN / 1000)<<" sum "<< (rv/1000) <<" rsum "<< (rsum/1000)<<std::endl;
        print_ll(ethis);
        print_ll(enext);
        std::cout<<std::endl;
        // distance left:
        edelta = eb - enext;
        double last_l = edelta.norm() * ELL_RMEAN;
        if (last_l < ELL_FLATSTEP)
        {
            std::cout<<"last step "<<(last_l/1000)<<std::endl;
            rv += last_l;
            break;
        }
    }
    return Py_BuildValue("d", rv);
}

static PyObject * // cartdist((na0,na1,na2), az, (nb0,nb1,nb2), bz) -> distance : return cartesian distance
ccnv_cartdist_ellipsoid(PyObject *self, PyObject *args)
{
    Eigen::Vector3d ea;
    Eigen::Vector3d eb;
    double a0, a1, a2, b0, b1, b2, az, bz;

    if (!PyArg_ParseTuple(args, "(ddd)d(ddd)d:cartdist_ellipse", &a0, &a1, &a2, &az, &b0, &b1, &b2, &bz))
        return NULL;

    ea[0] = a0; ea[1] = a1; ea[2] = a2;
    eb[0] = b0; eb[1] = b1; eb[2] = b2;

    Eigen::Vector3d ev = n_EA_E_and_n_EB_E2p_AB_E( ea, eb, az, bz, ELL_A, ELL_F );

    return Py_BuildValue("d", ev.norm());
}

static PyObject * // cartdist((na0,na1,na2), az, (nb0,nb1,nb2), bz) -> distance : return cartesian distance
ccnv_cartdist(PyObject *self, PyObject *args)
{
    Eigen::Vector3d ea;
    Eigen::Vector3d eb;
    double a0, a1, a2, b0, b1, b2, az, bz;

    if (!PyArg_ParseTuple(args, "(ddd)d(ddd)d:cartdist", &a0, &a1, &a2, &az, &b0, &b1, &b2, &bz))
        return NULL;

    ea[0] = a0; ea[1] = a1; ea[2] = a2;
    eb[0] = b0; eb[1] = b1; eb[2] = b2;

    Eigen::Vector3d ev = eb - ea ;

    return Py_BuildValue("d", ev.norm() * ELL_RMEAN);
}

static PyMethodDef ccnvMethods[] = {
    {"ellipsoid", ccnv_ellipsoid, METH_VARARGS,
        "ellipsoid(a, f)\n\n\
         Set ellipsoid paramerers (semimajor, flattening).\n\
         Default: WGS84, metres."},
    {"lat_long2n_E",  ccnv_lat_long2n_E, METH_VARARGS,
        "lat_long2n_E((lat,lon)) -> (n0,n1,n2)\n\n\
         Convert lan/lon to nvector."},
    {"latlon2n",  ccnv_lat_long2n_E, METH_VARARGS,
        "latlon2n((lat,lon)) -> (n0,n1,n2)\n\n\
         Convert lan/lon to nvector."},
    {"n_E2lat_long",  ccnv_n_E2lat_long, METH_VARARGS,
        "n_E2lat_long((n0,n1,n2)) -> (lat,lon)\n\n\
         Convert nvector to lan/lon."},
    {"n2latlon",  ccnv_n_E2lat_long, METH_VARARGS,
        "n2latlon((n0,n1,n2)) -> (lat,lon)\n\n\
        Convert nvector to lan/lon."},
    {"n_EA_E_and_n_EB_E2p_AB_E",  ccnv_n_EA_E_and_n_EB_E2p_AB_E, METH_VARARGS,
        "From two positions A and B, finds the delta vector."},
    {"ab2d",  ccnv_n_EA_E_and_n_EB_E2p_AB_E, METH_VARARGS,
        "ab2d((na0,na1,na2), (nb0,nb1,nb2), za, zb) -> (dx, dy, dz)\n\n\
         From two positions A and B, finds the delta vector."},
    {"n_EA_E_and_p_AB_E2n_EB_E",  ccnv_n_EA_E_and_p_AB_E2n_EB_E, METH_VARARGS,
        "From position A and delta, finds position B. (without B's depth)."},
    {"ad2b",  ccnv_n_EA_E_and_p_AB_E2n_EB_E, METH_VARARGS,
        "ab2d((na0,na1,na2), (px, py, pz), za) -> (nb0,nb1,nb2)\n\n\
         From position A and delta, finds position B. (without B's depth)."},
    {"n_EA_E_and_p_AB_E2n_EB_E_and_z_EB",  ccnv_n_EA_E_and_p_AB_E2n_EB_E_and_z_EB, METH_VARARGS,
        "From position A and delta, finds position B. (with B's depth)."},
    {"ad2bz",  ccnv_n_EA_E_and_p_AB_E2n_EB_E_and_z_EB, METH_VARARGS,
        "ab2d((na0,na1,na2), (px, py, pz), za) -> ((nb0,nb1,nb2), bz)\n\n\
         From position A and delta, finds position B. (with B's depth)."},
    {"sizevec", ccnv_sizevec, METH_VARARGS,
        "sizevec((x, y, z), l) ->(x, y, z)\n\n\
         Normalize vector, then multiply by l" },
    {"advance2d", ccnv_advance2d, METH_VARARGS,
        "advance2d((na0,na1,na2), (nb0,nb1,nb2), l) -> ((nc0, nc1, nc2), (lat, lon))\n\n\
         Advance na in direction of nb for l units over gc line." },
    {"gcdist", ccnv_gcdist, METH_VARARGS,
        "gcdist((na0,na1,na2), (nb0,nb1,nb2)) -> r\n\n\
         Return grand circle distance assuming sphere"},
    {"gcdist_ellipsoid", ccnv_gcdist_ellipsoid, METH_VARARGS,
        "gcdist_ellipsoid((na0,na1,na2), (nb0,nb1,nb2)) -> r\n\n\
         Return grand circle distance assuming ellipsoid (iterative)"},
    {"cartdist", ccnv_cartdist, METH_VARARGS,
        "cartdist((na0,na1,na2), (nb0,nb1,nb2)) -> r\n\n\
         Return cartesian distance assuming sphere"},
    {"cartdist_ellipsoid", ccnv_cartdist_ellipsoid, METH_VARARGS,
        "cartdist_ellipsoid((na0,na1,na2), (nb0,nb1,nb2)) -> r\n\n\
         Return cartesian distance assuming ellipsoid"},

    //{"",  ccnv_, METH_VARARGS, "."},

    {NULL, NULL, 0, NULL}        /* Sentinel */
};

static struct PyModuleDef ccnvmodule = {
   PyModuleDef_HEAD_INIT,
   "ccnvector",   /* name of module */
   NULL, /* module documentation, may be NULL */
   -1,       /* size of per-interpreter state of the module,
                or -1 if the module keeps state in global variables. */
   ccnvMethods
};

PyMODINIT_FUNC
PyInit_ccnvector(void)
{
    return PyModule_Create(&ccnvmodule);
}


