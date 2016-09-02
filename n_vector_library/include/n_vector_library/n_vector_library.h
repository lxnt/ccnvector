#ifndef __N_VECTOR_LIBRARY_H
#define __N_VECTOR_LIBRARY_H

#include <iostream>
#include <eigen3/Eigen/Dense>


namespace n_vector_library {


  /* Struct with position given as latitude and longitude.
   */
  struct LatLongPosition {
    double latitude;    // [rad], geodetic latitude
    double longitude;   // [rad]
  };


  /* Struct with position given as n-vector and depth.
   */
  struct NVectorPosition {
    Eigen::Vector3d n_EB_E;     // [rad], geodetic latitude
    double z_EB;                // [m]
  };


  /* Semi-major axis [meter], (equatorial radius) of the Earth.
   */
  const double WGS84_A = 6378137.0;

  /* Flattening of the Earth.
   */
  const double WGS84_F = 1/298.257223563;


  /* Converts latitude and longitude to n-vector.
   *
   * n-vector (n_E) is calculated from (geodetic) latitude and longitude.
   *
   * IN:
   * latitude:  [rad]     Geodetic latitude
   * longitude: [rad]
   *
   * OUT:
   * n_E:       [no unit] n-vector decomposed in E
   *
   * See also n_E2lat_long.
   *
   *
   * This file is part of NavLab and is available from www.navlab.net/nvector
   *
   * The content of this file is based on the following publication:
   *
   * Gade, K. (2010). A Nonsingular Horizontal Position Representation, The Journal
   * of Navigation, Volume 63, Issue 03, pp 395-417, July 2010.
   * (www.navlab.net/Publications/A_Nonsingular_Horizontal_Position_Representation.pdf)
   *
   * This paper should be cited in publications using this file.
   *
   * Copyright (c) 2015, Norwegian Defence Research Establishment (FFI)
   * All rights reserved.
   *
   * Redistribution and use in source and binary forms, with or without
   * modification, are permitted provided that the following conditions are met:
   *
   * 1. Redistributions of source code must retain the above publication
   * information, copyright notice, this list of conditions and the following disclaimer.
   *
   * 2. Redistributions in binary form must reproduce the above publication
   * information, copyright notice, this list of conditions and the following disclaimer
   * in the documentation and/or other materials provided with the distribution.
   *
   * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
   * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
   * PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS
   * BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
   * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
   * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
   * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
   * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
   * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
   * THE POSSIBILITY OF SUCH DAMAGE.
   *
   * Originated: 1999.02.23 Kenneth Gade, FFI
   * Modified:   2015.02.20 Kenneth Gade, FFI: Added possibility of using alternative E axes
   * Modified:   2016.02.09 Magnus Baksaas, FFI: Implemented the function i C++
   */
  Eigen::Vector3d lat_long2n_E( const double& latitude, const double& longitude );


  /* Converts n-vector to latitude and longitude.
   *
   * Geodetic latitude and longitude are calculated from n-vector (n_E).
   *
   * IN:
   * n_E:       [no unit] n-vector decomposed in E
   *
   * OUT:
   * latitude:  [rad]     Geodetic latitude
   * longitude: [rad]
   *
   * See also lat_long2n_E.
   *
   *
   * This file is part of NavLab and is available from www.navlab.net/nvector
   *
   * The content of this file is based on the following publication:
   *
   * Gade, K. (2010). A Nonsingular Horizontal Position Representation, The Journal
   * of Navigation, Volume 63, Issue 03, pp 395-417, July 2010.
   * (www.navlab.net/Publications/A_Nonsingular_Horizontal_Position_Representation.pdf)
   *
   * This paper should be cited in publications using this file.
   *
   * Copyright (c) 2015, Norwegian Defence Research Establishment (FFI)
   * All rights reserved.
   *
   * Redistribution and use in source and binary forms, with or without
   * modification, are permitted provided that the following conditions are met:
   *
   * 1. Redistributions of source code must retain the above publication
   * information, copyright notice, this list of conditions and the following disclaimer.
   *
   * 2. Redistributions in binary form must reproduce the above publication
   * information, copyright notice, this list of conditions and the following disclaimer
   * in the documentation and/or other materials provided with the distribution.
   *
   * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
   * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
   * PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS
   * BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
   * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
   * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
   * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
   * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
   * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
   * THE POSSIBILITY OF SUCH DAMAGE.
   *
   * Originated: 1999.02.23 Kenneth Gade, FFI
   * Modified:   2004.11.23 Kenneth Gade, FFI: Accepts vectorized input
   * Modified:   2015.02.20 Kenneth Gade, FFI: Added possibility of using alternative E axes
   * Modified:   2016.02.09 Magnus Baksaas, FFI: Implemented the function i C++
   */
  LatLongPosition n_E2lat_long( const Eigen::Vector3d& n_E );


  /* From two positions A and B, finds the delta position.
   *
   * The n-vectors for positions A (n_EA_E) and B (n_EB_E) are given. The
   * output is the delta vector from A to B (p_AB_E).
   * The calculation is exact, taking the ellipsity of the Earth into account.
   * It is also nonsingular as both n-vector and p-vector are nonsingular
   * (except for the center of the Earth).
   * The default ellipsoid model used is WGS-84, but other ellipsoids (or spheres)
   * might be specified.
   *
   * p_AB_E = n_EA_E_and_n_EB_E2p_AB_E(n_EA_E,n_EB_E,z_EA)
   * p_AB_E = n_EA_E_and_n_EB_E2p_AB_E(n_EA_E,n_EB_E,z_EA,z_EB)
   * Depth(s) of A, z_EA (and of B, z_EB) are also specified, z_EA = 0 (and z_EB = 0)
   * is used when not specified.
   *
   * p_AB_E = n_EA_E_and_n_EB_E2p_AB_E(n_EA_E,n_EB_E,z_EA,z_EB,a)
   * Spherical Earth with radius a is used instead of WGS-84.
   *
   * p_AB_E = n_EA_E_and_n_EB_E2p_AB_E(n_EA_E,n_EB_E,z_EA,z_EB,a,f)
   * Ellipsoidal Earth model with semi-major axis a and flattening f is used
   * instead of WGS-84.
   *
   * IN:
   * n_EA_E:  [no unit] n-vector of position A, decomposed in E.
   * n_EB_E:  [no unit] n-vector of position B, decomposed in E.
   * z_EA:    [m]       (Optional, assumed to be zero if not given) Depth of system A,
   *                    relative to the ellipsoid (z_EA = -height).
   * z_EB:    [m]       (Optional, assumed to be zero if not given) Depth of system B,
   *                    relative to the ellipsoid (z_EB = -height).
   * a:       [m]       (Optional) Semi-major axis of the Earth ellipsoid
   * f:       [no unit] (Optional) Flattening of the Earth ellipsoid
   *
   * OUT:
   * p_AB_E:  [m]       Position vector from A to B, decomposed in E.
   *
   * See also n_EA_E_and_p_AB_E2n_EB_E, p_EB_E2n_EB_E, n_EB_E2p_EB_E.
   *
   *
   * This file is part of NavLab and is available from www.navlab.net/nvector
   *
   * The content of this file is based on the following publication:
   *
   * Gade, K. (2010). A Nonsingular Horizontal Position Representation, The Journal
   * of Navigation, Volume 63, Issue 03, pp 395-417, July 2010.
   * (www.navlab.net/Publications/A_Nonsingular_Horizontal_Position_Representation.pdf)
   *
   * This paper should be cited in publications using this file.
   *
   * Copyright (c) 2015, Norwegian Defence Research Establishment (FFI)
   * All rights reserved.
   *
   * Redistribution and use in source and binary forms, with or without
   * modification, are permitted provided that the following conditions are met:
   *
   * 1. Redistributions of source code must retain the above publication
   * information, copyright notice, this list of conditions and the following disclaimer.
   *
   * 2. Redistributions in binary form must reproduce the above publication
   * information, copyright notice, this list of conditions and the following disclaimer
   * in the documentation and/or other materials provided with the distribution.
   *
   * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
   * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
   * PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS
   * BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
   * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
   * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
   * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
   * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
   * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
   * THE POSSIBILITY OF SUCH DAMAGE.
   *
   * Originated: 2004.07.07 Kenneth Gade, FFI
   * Modified:   2016.02.09 Magnus Baksaas, FFI: Implemented the function i C++
   */
  Eigen::Vector3d n_EA_E_and_n_EB_E2p_AB_E( const Eigen::Vector3d& n_EA_E,Eigen::Vector3d& n_EB_E, const double& z_EA = 0, const double& z_EB = 0, const double& a = WGS84_A, const double& f = WGS84_F );


  /* From position A and delta, finds position B.
   *
   * The n-vector for position A (n_EA_E) and the position-vector from position
   * A to position B (p_AB_E) are given. The output is the n-vector of position
   * B (n_EB_E) and depth of B (z_EB).
   * The calculation is exact, taking the ellipsity of the Earth into account.
   * It is also nonsingular as both n-vector and p-vector are nonsingular
   * (except for the center of the Earth).
   * The default ellipsoid model used is WGS-84, but other ellipsoids (or spheres)
   * might be specified.
   *
   * n_EB_E = n_EA_E_and_p_AB_E2n_EB_E(n_EA_E,p_AB_E,z_EA) Depth of A, z_EA,
   * is also specified, z_EA = 0 is used when not specified.
   *
   * n_EB_E = n_EA_E_and_p_AB_E2n_EB_E(n_EA_E,p_AB_E,z_EA,a)
   * Spherical Earth with radius a is used instead of WGS-84.
   *
   * n_EB_E = n_EA_E_and_p_AB_E2n_EB_E(n_EA_E,p_AB_E,z_EA,a,f)
   * Ellipsoidal Earth model with semi-major axis a and flattening f is used
   * instead of WGS-84.
   *
   * IN:
   * n_EA_E:  [no unit] n-vector of position A, decomposed in E.
   * p_AB_E:  [m]       Position vector from A to B, decomposed in E.
   * z_EA:    [m]       (Optional, assumed to be zero if not given) Depth of system A,
   *                    relative to the ellipsoid (z_EA = -height).
   * a:       [m]       (Optional) Semi-major axis of the Earth ellipsoid
   * f:       [no unit] (Optional) Flattening of the Earth ellipsoid
   *
   * OUT:
   * n_EB_E:  [no unit] n-vector of position B, decomposed in E.
   *
   * See also n_EA_E_and_n_EB_E2p_AB_E, p_EB_E2n_EB_E, n_EB_E2p_EB_E.
   *
   *
   * This file is part of NavLab and is available from www.navlab.net/nvector
   *
   * The content of this file is based on the following publication:
   *
   * Gade, K. (2010). A Nonsingular Horizontal Position Representation, The Journal
   * of Navigation, Volume 63, Issue 03, pp 395-417, July 2010.
   * (www.navlab.net/Publications/A_Nonsingular_Horizontal_Position_Representation.pdf)
   *
   * This paper should be cited in publications using this file.
   *
   * Copyright (c) 2015, Norwegian Defence Research Establishment (FFI)
   * All rights reserved.
   *
   * Redistribution and use in source and binary forms, with or without
   * modification, are permitted provided that the following conditions are met:
   *
   * 1. Redistributions of source code must retain the above publication
   * information, copyright notice, this list of conditions and the following disclaimer.
   *
   * 2. Redistributions in binary form must reproduce the above publication
   * information, copyright notice, this list of conditions and the following disclaimer
   * in the documentation and/or other materials provided with the distribution.
   *
   * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
   * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
   * PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS
   * BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
   * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
   * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
   * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
   * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
   * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
   * THE POSSIBILITY OF SUCH DAMAGE.
   *
   * Originated: 2004.07.07 Kenneth Gade, FFI
   * Modified:   2016.02.09 Magnus Baksaas, FFI: Implemented the function i C++
   */
  Eigen::Vector3d n_EA_E_and_p_AB_E2n_EB_E( const Eigen::Vector3d& n_EA_E, const Eigen::Vector3d& p_AB_E, const double& z_EA = 0, const double& a = WGS84_A, const double& f = WGS84_F );


  /* From position A and delta, finds position B.
   *
   * The n-vector for position A (n_EA_E) and the position-vector from position
   * A to position B (p_AB_E) are given. The output is the n-vector of position
   * B (n_EB_E) and depth of B (z_EB).
   * The calculation is exact, taking the ellipsity of the Earth into account.
   * It is also nonsingular as both n-vector and p-vector are nonsingular
   * (except for the center of the Earth).
   * The default ellipsoid model used is WGS-84, but other ellipsoids (or spheres)
   * might be specified.
   *
   * [n_EB_E,z_EB] = n_EA_E_and_p_AB_E2n_EB_E(n_EA_E,p_AB_E,z_EA) Depth of A, z_EA,
   * is also specified, z_EA = 0 is used when not specified.
   *
   * [n_EB_E,z_EB] = n_EA_E_and_p_AB_E2n_EB_E(n_EA_E,p_AB_E,z_EA,a)
   * Spherical Earth with radius a is used instead of WGS-84.
   *
   * [n_EB_E,z_EB] = n_EA_E_and_p_AB_E2n_EB_E(n_EA_E,p_AB_E,z_EA,a,f)
   * Ellipsoidal Earth model with semi-major axis a and flattening f is used
   * instead of WGS-84.
   *
   * IN:
   * n_EA_E:  [no unit] n-vector of position A, decomposed in E.
   * p_AB_E:  [m]       Position vector from A to B, decomposed in E.
   * z_EA:    [m]       (Optional, assumed to be zero if not given) Depth of system A,
   *                    relative to the ellipsoid (z_EA = -height).
   * a:       [m]       (Optional) Semi-major axis of the Earth ellipsoid
   * f:       [no unit] (Optional) Flattening of the Earth ellipsoid
   *
   * OUT:
   * n_EB_E:  [no unit] n-vector of position B, decomposed in E.
   * z_EB:    [m]       Depth of system B, relative to the ellipsoid (z_EB = -height).
   *
   * See also n_EA_E_and_n_EB_E2p_AB_E, p_EB_E2n_EB_E, n_EB_E2p_EB_E.
   *
   *
   * This file is part of NavLab and is available from www.navlab.net/nvector
   *
   * The content of this file is based on the following publication:
   *
   * Gade, K. (2010). A Nonsingular Horizontal Position Representation, The Journal
   * of Navigation, Volume 63, Issue 03, pp 395-417, July 2010.
   * (www.navlab.net/Publications/A_Nonsingular_Horizontal_Position_Representation.pdf)
   *
   * This paper should be cited in publications using this file.
   *
   * Copyright (c) 2015, Norwegian Defence Research Establishment (FFI)
   * All rights reserved.
   *
   * Redistribution and use in source and binary forms, with or without
   * modification, are permitted provided that the following conditions are met:
   *
   * 1. Redistributions of source code must retain the above publication
   * information, copyright notice, this list of conditions and the following disclaimer.
   *
   * 2. Redistributions in binary form must reproduce the above publication
   * information, copyright notice, this list of conditions and the following disclaimer
   * in the documentation and/or other materials provided with the distribution.
   *
   * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
   * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
   * PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS
   * BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
   * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
   * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
   * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
   * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
   * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
   * THE POSSIBILITY OF SUCH DAMAGE.
   *
   * Originated: 2004.07.07 Kenneth Gade, FFI
   * Modified:   2016.02.09 Magnus Baksaas, FFI: Implemented the function i C++
   */
  NVectorPosition n_EA_E_and_p_AB_E2n_EB_E_and_z_EB( const Eigen::Vector3d& n_EA_E, const Eigen::Vector3d& p_AB_E, const double& z_EA = 0, const double& a = WGS84_A, const double& f = WGS84_F );


  /* Converts n-vector to Cartesian position vector in meters.
   *
   * The position of B (typically body) relative to E (typically Earth) is
   * given into this function as n-vector, n_EB_E. The function converts
   * to cartesian position vector ("ECEF-vector"), p_EB_E, in meters.
   * The calculation is exact, taking the ellipsity of the Earth into account.
   * It is also nonsingular as both n-vector and p-vector are nonsingular
   * (except for the center of the Earth).
   * The default ellipsoid model used is WGS-84, but other ellipsoids (or spheres)
   * might be specified.
   *
   * p_EB_E = n_EB_E2p_EB_E(n_EB_E,z_EB) Depth of B, z_EB, is also specified,
   * z_EB = 0 is used when not specified.
   *
   * p_EB_E = n_EB_E2p_EB_E(n_EB_E,z_EB,a) Spherical Earth with radius a is
   * used instead of WGS-84.
   *
   * p_EB_E = n_EB_E2p_EB_E(n_EB_E,z_EB,a,f) Ellipsoidal Earth model with
   * semi-major axis a and flattening f is used instead of WGS-84.
   *
   *
   * IN:
   * n_EB_E:  [no unit] n-vector of position B, decomposed in E.
   * z_EB:    [m]       (Optional, assumed to be zero if not given) Depth of system B,
   *                    relative to the ellipsoid (z_EB = -height)
   * a:       [m]       (Optional) Semi-major axis of the Earth ellipsoid
   * f:       [no unit] (Optional) Flattening of the Earth ellipsoid
   *
   * OUT:
   * p_EB_E:  [m]       Cartesian position vector from E to B, decomposed in E.
   *
   * See also p_EB_E2n_EB_E, n_EA_E_and_p_AB_E2n_EB_E, n_EA_E_and_n_EB_E2p_AB_E.
   *
   *
   * This file is part of NavLab and is available from www.navlab.net/nvector
   *
   * The content of this file is based on the following publication:
   *
   * Gade, K. (2010). A Nonsingular Horizontal Position Representation, The Journal
   * of Navigation, Volume 63, Issue 03, pp 395-417, July 2010.
   * (www.navlab.net/Publications/A_Nonsingular_Horizontal_Position_Representation.pdf)
   *
   * This paper should be cited in publications using this file.
   *
   * Copyright (c) 2015, Norwegian Defence Research Establishment (FFI)
   * All rights reserved.
   *
   * Redistribution and use in source and binary forms, with or without
   * modification, are permitted provided that the following conditions are met:
   *
   * 1. Redistributions of source code must retain the above publication
   * information, copyright notice, this list of conditions and the following disclaimer.
   *
   * 2. Redistributions in binary form must reproduce the above publication
   * information, copyright notice, this list of conditions and the following disclaimer
   * in the documentation and/or other materials provided with the distribution.
   *
   * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
   * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
   * PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS
   * BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
   * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
   * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
   * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
   * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
   * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
   * THE POSSIBILITY OF SUCH DAMAGE.
   *
   * Originated: 2004.11.17 Kenneth Gade and Brita Hafskjold, FFI
   * Modified:   2015.02.20 Kenneth Gade, FFI: Added possibility of using alternative E axes
   * Modified:   2016.02.09 Magnus Baksaas, FFI: Implemented the function i C++
   */
  Eigen::Vector3d n_EB_E2p_EB_E( const Eigen::Vector3d& n_EB_E, const double& z_EB = 0, const double& a = WGS84_A, const double& f = WGS84_F );


  /* Converts Cartesian position vector in meters to n-vector.
   *
   * The position of B (typically body) relative to E (typically Earth) is
   * given into this function as cartesian position vector p_EB_E, in meters
   * ("ECEF-vector"). The function converts to n-vector, n_EB_E and its
   * depth, z_EB.
   * The calculation is exact, taking the ellipsity of the Earth into account.
   * It is also nonsingular as both n-vector and p-vector are nonsingular
   * (except for the center of the Earth).
   * The default ellipsoid model used is WGS-84, but other ellipsoids (or spheres)
   * might be specified.
   *
   * n_EB_E = p_EB_E2n_EB_E(p_EB_E,a) Spherical Earth with radius a is
   * used instead of WGS-84.
   *
   * n_EB_E = p_EB_E2n_EB_E(p_EB_E,a,f) Ellipsoidal Earth model with
   * semi-major axis a and flattening f is used instead of WGS-84.
   *
   *
   * IN:
   * p_EB_E: [m]       Cartesian position vector from E to B, decomposed in E.
   * a:      [m]       (Optional) Semi-major axis of the Earth ellipsoid
   * f:      [no unit] (Optional) Flattening of the Earth ellipsoid
   *
   * OUT:
   * n_EB_E: [no unit] n-vector  representation of position B, decomposed in E.
   * z_EB:   [m]       Depth of system B relative to the ellipsoid (z_EB = -height).
   *
   * See also n_EB_E2p_EB_E, n_EA_E_and_p_AB_E2n_EB_E, n_EA_E_and_n_EB_E2p_AB_E.
   *
   *
   * This file is part of NavLab and is available from www.navlab.net/nvector
   *
   * The content of this file is based on the following publication:
   *
   * Gade, K. (2010). A Nonsingular Horizontal Position Representation, The Journal
   * of Navigation, Volume 63, Issue 03, pp 395-417, July 2010.
   * (www.navlab.net/Publications/A_Nonsingular_Horizontal_Position_Representation.pdf)
   *
   * This paper should be cited in publications using this file.
   *
   * Copyright (c) 2015, Norwegian Defence Research Establishment (FFI)
   * All rights reserved.
   *
   * Redistribution and use in source and binary forms, with or without
   * modification, are permitted provided that the following conditions are met:
   *
   * 1. Redistributions of source code must retain the above publication
   * information, copyright notice, this list of conditions and the following disclaimer.
   *
   * 2. Redistributions in binary form must reproduce the above publication
   * information, copyright notice, this list of conditions and the following disclaimer
   * in the documentation and/or other materials provided with the distribution.
   *
   * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
   * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
   * PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS
   * BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
   * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
   * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
   * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
   * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
   * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
   * THE POSSIBILITY OF SUCH DAMAGE.
   *
   * Originated: 2004.11.17 Kenneth Gade and Brita Hafskjold, FFI
   * Modified:   2007.03.02 Brita Hafskjold Gade, FFI:
   *               Replaced formulas to get full numerical accuracy at all positions
   * Modified:   2014.08.22 Kenneth Gade, FFI:
   *               Added possibility of vectorized input/output
   * Modified:   2016.02.09 Magnus Baksaas, FFI: Implemented the function i C++
   */
  Eigen::Vector3d p_EB_E2n_EB_E( const Eigen::Vector3d& p_EB_E, const double& a = WGS84_A, const double& f = WGS84_F );


  /* Converts Cartesian position vector in meters to n-vector.
   *
   * The position of B (typically body) relative to E (typically Earth) is
   * given into this function as cartesian position vector p_EB_E, in meters
   * ("ECEF-vector"). The function converts to n-vector, n_EB_E and its
   * depth, z_EB.
   * The calculation is exact, taking the ellipsity of the Earth into account.
   * It is also nonsingular as both n-vector and p-vector are nonsingular
   * (except for the center of the Earth).
   * The default ellipsoid model used is WGS-84, but other ellipsoids (or spheres)
   * might be specified.
   *
   * [n_EB_E,z_EB] = p_EB_E2n_EB_E(p_EB_E,a) Spherical Earth with radius a is
   * used instead of WGS-84.
   *
   * [n_EB_E,z_EB] = p_EB_E2n_EB_E(p_EB_E,a,f) Ellipsoidal Earth model with
   * semi-major axis a and flattening f is used instead of WGS-84.
   *
   *
   * IN:
   * p_EB_E: [m]       Cartesian position vector from E to B, decomposed in E.
   * a:      [m]       (Optional) Semi-major axis of the Earth ellipsoid
   * f:      [no unit] (Optional) Flattening of the Earth ellipsoid
   *
   * OUT:
   * n_EB_E: [no unit] n-vector  representation of position B, decomposed in E.
   * z_EB:   [m]       Depth of system B relative to the ellipsoid (z_EB = -height).
   *
   * See also n_EB_E2p_EB_E, n_EA_E_and_p_AB_E2n_EB_E, n_EA_E_and_n_EB_E2p_AB_E.
   *
   *
   * This file is part of NavLab and is available from www.navlab.net/nvector
   *
   * The content of this file is based on the following publication:
   *
   * Gade, K. (2010). A Nonsingular Horizontal Position Representation, The Journal
   * of Navigation, Volume 63, Issue 03, pp 395-417, July 2010.
   * (www.navlab.net/Publications/A_Nonsingular_Horizontal_Position_Representation.pdf)
   *
   * This paper should be cited in publications using this file.
   *
   * Copyright (c) 2015, Norwegian Defence Research Establishment (FFI)
   * All rights reserved.
   *
   * Redistribution and use in source and binary forms, with or without
   * modification, are permitted provided that the following conditions are met:
   *
   * 1. Redistributions of source code must retain the above publication
   * information, copyright notice, this list of conditions and the following disclaimer.
   *
   * 2. Redistributions in binary form must reproduce the above publication
   * information, copyright notice, this list of conditions and the following disclaimer
   * in the documentation and/or other materials provided with the distribution.
   *
   * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
   * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
   * PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS
   * BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
   * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
   * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
   * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
   * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
   * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
   * THE POSSIBILITY OF SUCH DAMAGE.
   *
   * Originated: 2004.11.17 Kenneth Gade and Brita Hafskjold, FFI
   * Modified:   2007.03.02 Brita Hafskjold Gade, FFI:
   *               Replaced formulas to get full numerical accuracy at all positions
   * Modified:   2014.08.22 Kenneth Gade, FFI:
   *               Added possibility of vectorized input/output
   * Modified:   2016.02.09 Magnus Baksaas, FFI: Implemented the function i C++
   */
  NVectorPosition p_EB_E2n_EB_E_and_z_EB( const Eigen::Vector3d& p_EB_E, const double& a = WGS84_A, const double& f = WGS84_F );


  /* Finds n-vector from R_EN.
   *
   * n-vector is found from the rotation matrix (direction cosine matrix)
   * R_EN.
   *
   * IN:
   * R_EN:  [no unit] Rotation matrix (direction cosine matrix)
   *
   * OUT:
   * n_E:   [no unit] n-vector decomposed in E
   *
   * See also n_E2R_EN.
   *
   *
   * This file is part of NavLab and is available from www.navlab.net/nvector
   *
   * The content of this file is based on the following publication:
   *
   * Gade, K. (2010). A Nonsingular Horizontal Position Representation, The Journal
   * of Navigation, Volume 63, Issue 03, pp 395-417, July 2010.
   * (www.navlab.net/Publications/A_Nonsingular_Horizontal_Position_Representation.pdf)
   *
   * This paper should be cited in publications using this file.
   *
   * Copyright (c) 2015, Norwegian Defence Research Establishment (FFI)
   * All rights reserved.
   *
   * Redistribution and use in source and binary forms, with or without
   * modification, are permitted provided that the following conditions are met:
   *
   * 1. Redistributions of source code must retain the above publication
   * information, copyright notice, this list of conditions and the following disclaimer.
   *
   * 2. Redistributions in binary form must reproduce the above publication
   * information, copyright notice, this list of conditions and the following disclaimer
   * in the documentation and/or other materials provided with the distribution.
   *
   * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
   * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
   * PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS
   * BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
   * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
   * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
   * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
   * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
   * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
   * THE POSSIBILITY OF SUCH DAMAGE.
   *
   * Originated: 1999.02.23 Kenneth Gade, FFI
   * Modified:   2016.02.09 Magnus Baksaas, FFI: Implemented the function i C++
   */
  Eigen::Vector3d R_EN2n_E( const Eigen::Matrix3d& R_EN );


  /* Finds the rotation matrix R_EN from n-vector.
   *
   * The rotation matrix (direction cosine matrix) R_EN is calculated based
   * on n-vector (n_E).
   *
   * IN:
   * n_E:   [no unit] n-vector decomposed in E
   *
   * OUT:
   * R_EN:  [no unit] The resulting rotation matrix (direction cosine matrix)
   *
   * See also R_EN2n_E.
   *
   *
   * This file is part of NavLab and is available from www.navlab.net/nvector
   *
   * The content of this file is based on the following publication:
   *
   * Gade, K. (2010). A Nonsingular Horizontal Position Representation, The Journal
   * of Navigation, Volume 63, Issue 03, pp 395-417, July 2010.
   * (www.navlab.net/Publications/A_Nonsingular_Horizontal_Position_Representation.pdf)
   *
   * This paper should be cited in publications using this file.
   *
   * Copyright (c) 2015, Norwegian Defence Research Establishment (FFI)
   * All rights reserved.
   *
   * Redistribution and use in source and binary forms, with or without
   * modification, are permitted provided that the following conditions are met:
   *
   * 1. Redistributions of source code must retain the above publication
   * information, copyright notice, this list of conditions and the following disclaimer.
   *
   * 2. Redistributions in binary form must reproduce the above publication
   * information, copyright notice, this list of conditions and the following disclaimer
   * in the documentation and/or other materials provided with the distribution.
   *
   * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
   * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
   * PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS
   * BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
   * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
   * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
   * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
   * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
   * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
   * THE POSSIBILITY OF SUCH DAMAGE.
   *
   * Originated: 2015.02.23 Kenneth Gade, FFI
   * Modified:   2016.02.09 Magnus Baksaas, FFI: Implemented the function i C++
   */
  Eigen::Matrix3d n_E2R_EN( const Eigen::Vector3d& n_E );


}

#endif