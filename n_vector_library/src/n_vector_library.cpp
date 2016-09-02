#include <iostream>
#include <eigen3/Eigen/Dense>
#include "n_vector_library/n_vector_library.h"


namespace n_vector_library {


  Eigen::Vector3d lat_long2n_E( const double& latitude, const double& longitude ) {

    double sin_lat = sin(latitude);
    double cos_lat = cos(latitude);
    double sin_lon = sin(longitude);
    double cos_lon = cos(longitude);

    Eigen::Vector3d n_E;
    n_E << cos_lat*cos_lon,
    cos_lat*sin_lon,
    sin_lat;

    return n_E;
  }


  LatLongPosition n_E2lat_long( const Eigen::Vector3d& n_E ) {

    LatLongPosition position;
    position.latitude  = atan2( n_E(2), sqrt( n_E(0)*n_E(0) + n_E(1)*n_E(1) ) );
    position.longitude = atan2( n_E(1), n_E(0) );

    return position;
  }


  Eigen::Vector3d n_EA_E_and_n_EB_E2p_AB_E( const Eigen::Vector3d& n_EA_E,Eigen::Vector3d& n_EB_E, const double& z_EA, const double& z_EB, const double& a, const double& f ) {

    Eigen::Vector3d p_EA_E = n_EB_E2p_EB_E( n_EA_E, z_EA, a, f );
    Eigen::Vector3d p_EB_E = n_EB_E2p_EB_E( n_EB_E, z_EB, a, f );

    return p_EB_E - p_EA_E;
  }


  Eigen::Vector3d n_EA_E_and_p_AB_E2n_EB_E( const Eigen::Vector3d& n_EA_E, const Eigen::Vector3d& p_AB_E, const double& z_EA, const double& a, const double& f ) {

    Eigen::Vector3d p_EA_E = n_EB_E2p_EB_E( n_EA_E, z_EA, a, f);
    Eigen::Vector3d p_EB_E = p_EA_E + p_AB_E;

    return p_EB_E2n_EB_E( p_EB_E, a, f);
  }


  NVectorPosition n_EA_E_and_p_AB_E2n_EB_E_and_z_EB( const Eigen::Vector3d& n_EA_E, const Eigen::Vector3d& p_AB_E, const double& z_EA, const double& a, const double& f ) {

    Eigen::Vector3d p_EA_E = n_EB_E2p_EB_E( n_EA_E, z_EA, a, f);
    Eigen::Vector3d p_EB_E = p_EA_E + p_AB_E;

    return p_EB_E2n_EB_E_and_z_EB( p_EB_E, a, f);
  }


  Eigen::Vector3d n_EB_E2p_EB_E( const Eigen::Vector3d& n_EB_E, const double& z_EB, const double& a, const double& f ) {

    // Ensures unit length
    Eigen::Vector3d nn_EB_E = n_EB_E;
    nn_EB_E.normalize();

    // semi-minor axis:
    double b = a * (1 - f);

    // The following code implements equation (22) in Gade (2010):
    double denominator = sqrt( a*a/(b*b)*nn_EB_E(0)*nn_EB_E(0) + a*a/(b*b)*nn_EB_E(1)*nn_EB_E(1) + nn_EB_E(2)*nn_EB_E(2) );

    // We first calculate the position at the origin of coordinate system L,
    // which has the same n-vector as B (n_EL_E = n_EB_E),
    // but lies at the surface of the Earth (z_EL = 0).
    Eigen::Vector3d p_EB_E;
    p_EB_E << a*a/(b*b)*nn_EB_E(0),
    a*a/(b*b)*nn_EB_E(1),
    nn_EB_E(2);
    p_EB_E *= b/denominator;

    p_EB_E -= z_EB*n_EB_E;
    return p_EB_E;
  }


  Eigen::Vector3d p_EB_E2n_EB_E( const Eigen::Vector3d& p_EB_E, const double& a, const double& f )  {

    double a2 = a*a;
    double e2 = 2*f - f*f;
    double e4 = e2*e2;

    double q = ( ( 1 - e2 )/a2 )*p_EB_E(2)*p_EB_E(2);
    double p = ( p_EB_E(0)*p_EB_E(0) + p_EB_E(1)*p_EB_E(1) )/a2;
    double r = ( p + q - e4 )/6;
    double s = e4*p*q/(4*r*r*r);
    double t = cbrt( 1 + s + sqrt( s*( 2 + s ) ) );
    double u = r*( 1 + t + 1/t );
    double v = sqrt( u*u + e4*q );
    double w = e2*( u + v - q )/(2*v);
    double k = sqrt( u + v + w*w ) - w;
    double d = k*sqrt( p_EB_E(0)*p_EB_E(0) + p_EB_E(1)*p_EB_E(1) )/( k + e2 );

    double scalar = 1/( sqrt( d*d + p_EB_E(2)*p_EB_E(2) ) );

    Eigen::Vector3d n_E;
    n_E << scalar*k/( k + e2 )*p_EB_E(0),
    scalar*k/( k + e2 )*p_EB_E(1),
    scalar*p_EB_E(2);

    n_E.normalize();

    return n_E;
  }


  NVectorPosition p_EB_E2n_EB_E_and_z_EB( const Eigen::Vector3d& p_EB_E, const double& a, const double& f ) {

    double a2 = a*a;
    double e2 = 2*f - f*f;
    double e4 = e2*e2;

    double q = ( ( 1 - e2 )/a2 )*p_EB_E(2)*p_EB_E(2);
    double p = ( p_EB_E(0)*p_EB_E(0) + p_EB_E(1)*p_EB_E(1) )/a2;
    double r = ( p + q - e4 )/6;
    double s = e4*p*q/(4*r*r*r);
    double t = cbrt( 1 + s + sqrt( s*( 2 + s ) ) );
    double u = r*( 1 + t + 1/t );
    double v = sqrt( u*u + e4*q );
    double w = e2*( u + v - q )/(2*v);
    double k = sqrt( u + v + w*w ) - w;
    double d = k*sqrt( p_EB_E(0)*p_EB_E(0) + p_EB_E(1)*p_EB_E(1) )/( k + e2 );

    double scalar = 1/( sqrt( d*d + p_EB_E(2)*p_EB_E(2) ) );

    Eigen::Vector3d n_E;
    n_E << scalar*k/( k + e2 )*p_EB_E(0),
    scalar*k/( k + e2 )*p_EB_E(1),
    scalar*p_EB_E(2);

    n_E.normalize();

    NVectorPosition position;
    position.n_EB_E = n_E;
    position.z_EB   = -( k + e2 - 1 )/k * sqrt( d*d + p_EB_E(2)*p_EB_E(2) );

    return position;
  }


  Eigen::Vector3d R_EN2n_E( const Eigen::Matrix3d& R_EN ) {
    Eigen::Vector3d unit( 0, 0, -1 );
    return R_EN*unit;
  };


  Eigen::Matrix3d n_E2R_EN( const Eigen::Vector3d& n_E ) {

    Eigen::Vector3d z_E(0,0,1);

    Eigen::Vector3d Nx_E;
    Eigen::Vector3d Ny_E;
    Eigen::Vector3d Ny_E_direction;
    Eigen::Vector3d Nz_E;

    Eigen::Matrix3d R_EN;

    // Ensures unit length
    Eigen::Vector3d nn_E = n_E;
    nn_E.normalize();

    // Find z-axis of N (Nz):
    Nz_E = -nn_E;                         // z-axis of N (down) points opposite to n-vector

    // Find y-axis of N (East) (remember that N is singular at Poles)
    Ny_E_direction = z_E.cross(nn_E);     // Ny points perpendicular to the plane formed by n-vector and Earth's spin axis

    if( Ny_E_direction.norm() != 0 ) {    //  outside Poles
      Ny_E_direction.normalize();
      Ny_E = Ny_E_direction;
    }
    else {                                // Pole position:
      Ny_E << 0, 1, 0;
    }

    // Find x-axis of N (North):
    Nx_E = Ny_E.cross(Nz_E);              // Final axis found by right hand rule

    R_EN.col(0) = Nx_E;
    R_EN.col(1) = Ny_E;
    R_EN.col(2) = Nz_E;

    return R_EN;
  }


}