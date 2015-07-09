/* $URL: https://repos.ser.asu.edu/svn/lroc/ruby_spice/tags/release_20140107.1/spice.c $
*
* Copyright (C) 2013 by Arizona State University and Mark Robinson.
* All rights reserved.
*
* This file has been released under the modified BSD license.
* See COPYING in the distribution package for details.
*
* Author: Nick Estes <nme@ser.asu.edu>
*
* $Author: alicht $
*   $Date: 2013-10-25 13:24:11 -0700 (Fri, 25 Oct 2013) $
*    $Rev: 12925 $
*/

#include "ruby.h"
#include <stdbool.h>
#include "SpiceUsr.h"
#include "signal.h"

/** 
 * \mainpage Ruby Spice
 * This is a partial wrapper of the CSPICE library provided by 
 * <a href="http://naif.jpl.nasa.gov/naif/">NAIF</a>
 * It is not yet feature complete, but most of the useful functions
 * are present.
 *
 * This gem has been released under the modified BSD license.
 * See COPYING in the distribution package for details.
 *
 * Copyright (c) 2009-2013, Arizona State University and Mark Robinson
 * All rights reserved.
 *
 */

/** @file spice.c
 * This is the heart of ruby_spice.
 * All of the function wrappers tale place here.
 *
 */

/*Since doxygen can't see ruby stuff we'll include the skh here*/
/*!
 * @page "Spice Kernel Helper"
 * This class is used to parse out mk files
 * LROC mk files as a whole wayyy blow past the 1000 file handle limit
 * This pares it down to only the kernels which we need.
 * This file is already 'required' when you run 'require "spice"'
 *
 *
 * Example Usage:
 *@code
 * skh = SpiceKernelHelper.new
 *
 * #Loop through an array of kernels, furnshing
 * #non-mk kernels directly
 * kernels.each do |kernel|
 *   if kernel =~ /mk$/
 *     skh.parse kernel
 *   else
 *     furnsh kernel
 *   end
 * end
 * furnsh *skh.base_kernels
 * furnsh *skh.detail_kernels(et2rb(start_et))
 * @endcode
 *
 */

/*!
 *@page "A Note Regarding Aberration Corrections"
 *abcorr      indicates the aberration corrections to be applied when
               computing the target's position and orientation.
 
               For remote sensing applications, where the apparent 
               sub-observer point seen by the observer is desired, 
               normally either of the corrections  
             
                  "LT+S"  
                  "CN+S" 
    
               should be used. These and the other supported options 
               are described below. `abcorr' may be any of the  
               following: 
 
                  "NONE"     Apply no correction. Return the  
                             geometric sub-observer point on the 
                             target body. 
 
               Let `lt' represent the one-way light time between the 
               observer and the sub-observer point (note: NOT 
               between the observer and the target body's center). 
               The following values of `abcorr' apply to the 
               "reception" case in which photons depart from the 
               sub-observer point's location at the light-time 
               corrected epoch et-lt and *arrive* at the observer's 
               location at `et': 
 
 
                  "LT"       Correct for one-way light time (also 
                             called "planetary aberration") using a 
                             Newtonian formulation. This correction 
                             yields the location of sub-observer 
                             point at the moment it emitted photons 
                             arriving at the observer at `et'. 
  
                             The light time correction uses an 
                             iterative solution of the light time 
                             equation. The solution invoked by the 
                             "LT" option uses one iteration. 
 
                             Both the target position as seen by the 
                             observer, and rotation of the target 
                             body, are corrected for light time. 
 
                  "LT+S"     Correct for one-way light time and stellar
                             aberration using a Newtonian formulation.
                             This option modifies the sub-observer
                             point obtained with the "LT" option to
                             account for the observer's velocity
                             relative to the solar system barycenter.
                             These corrections yield the apparent
                             sub-observer point.

                  "CN"       Converged Newtonian light time 
                             correction. In solving the light time 
                             equation, the "CN" correction iterates 
                             until the solution converges. Both the 
                             position and rotation of the target 
                             body are corrected for light time. 
 
                  "CN+S"     Converged Newtonian light time and 
                             stellar aberration corrections. This 
                             option produces a solution that is at 
                             least as accurate at that obtainable 
                             with the "LT+S" option. Whether the "CN+S" 
                             solution is substantially more accurate 
                             depends on the geometry of the 
                             participating objects and on the 
                             accuracy of the input data. In all 
                             cases this routine will execute more 
                             slowly when a converged solution is 
                             computed. 
                              
 
               The following values of `abcorr' apply to the 
               "transmission" case in which photons *depart* from 
               the observer's location at `et' and arrive at the 
               sub-observer point at the light-time corrected epoch 
               et+lt: 
 
                  "XLT"      "Transmission" case: correct for 
                             one-way light time using a Newtonian 
                             formulation. This correction yields the 
                             sub-observer location at the moment it 
                             receives photons emitted from the 
                             observer's location at `et'.  
 
                             The light time correction uses an 
                             iterative solution of the light time 
                             equation. The solution invoked by the 
                             "LT" option uses one iteration. 
 
                             Both the target position as seen by the 
                             observer, and rotation of the target 
                             body, are corrected for light time. 
 
                  "XLT+S"    "Transmission" case: correct for 
                             one-way light time and stellar 
                             aberration using a Newtonian 
                             formulation  This option modifies the 
                             sub-observer point obtained with the 
                             "XLT" option to account for the 
                             observer's velocity relative to the 
                             solar system barycenter. 
 
                  "XCN"      Converged Newtonian light time 
                             correction. This is the same as "XLT"
                             correction but with further iterations 
                             to a converged Newtonian light time 
                             solution.  
 
                  "XCN+S"    "Transmission" case: converged  
                             Newtonian light time and stellar  
                             aberration corrections. 
 
                Neither case nor white space are significant in
                `abcorr'. For example, the string

                  'Lt + s'

                is valid.

 *
 */

VALUE rb_eSpiceError;

sigset_t old_signal_mask;

void block_signals() {
  sigset_t mask_all;
  sigfillset(&mask_all);
  sigprocmask(SIG_BLOCK, &mask_all, &old_signal_mask);
}

void restore_signals() {
  sigprocmask(SIG_SETMASK, &old_signal_mask, NULL);
}

void check_spice_error() {
  SpiceInt lenout = 1024;
  char mShort[lenout], mExplain[lenout], mLong[lenout];

  if (failed_c()) {
    getmsg_c("SHORT", lenout, mShort);
    getmsg_c("EXPLAIN", lenout, mExplain);
    getmsg_c("LONG", lenout, mLong);
    reset_c();
    rb_raise(rb_eSpiceError, "%s\n%s\n\n%s\n", mShort, mExplain, mLong);
  }
}
/*!@fn VALUE furnsh(int argc, VALUE *argv, VALUE self)
 * @brief Load one or more SPICE kernels into a program.
 * @brief http://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/furnsh_c.html
 * @param Path to a Spice Kernel, can be a single kernel of .mk file
 *
 *@return returns a fixnum of the number of kernels loaded
 */
VALUE furnsh(int argc, VALUE *argv, VALUE self) {
  SpiceInt i;
  SpiceInt numkernels;

  block_signals();

  if (argc == 0) {
    rb_raise(rb_eArgError, "furnsh needs kernels!");
  } else {
    for(i=0; i < argc; i++) {
      Check_Type(argv[i], T_STRING);
    }
  }

  for(i=0; i < argc; i++) {
    furnsh_c(StringValuePtr(argv[i]));
  }

  ktotal_c("ALL", &numkernels);

  restore_signals();

  check_spice_error();

  return INT2FIX(numkernels);
}
/*! @fn VALUE unload(int argc, VALUE *argv, VALUE self)
 *@brief  Unload a SPICE kernel
 *@brief http://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/unload_c.html
 * @param Path to a Spice Kernel, can be a single kernel of .mk file
 *  *@return returns a fixnum of the number of kernels unloaded
 *
 */
VALUE unload(int argc, VALUE *argv, VALUE self) {
  SpiceInt i;
  SpiceInt numkernels;

  block_signals();

  if (argc == 0) {
    rb_raise(rb_eArgError, "unload needs kernels!");
  } else {
    for(i=0; i < argc; i++) {
      Check_Type(argv[i], T_STRING);
    }
  }

  for(i=0; i < argc; i++) {
    unload_c(StringValuePtr(argv[i]));
  }

  ktotal_c("ALL", &numkernels);

  restore_signals();

  check_spice_error();

  return INT2FIX(numkernels);
}

/*!@fn VALUE kclear(VALUE self)
 * @brief Clear the KEEPER system:  unload all kernels, clear the kernel pool, and re-initialize the system. 
 * @brief http://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/kclear_c.html
 * @param NONE
 *
 * @return true if sucessful
 */
VALUE kclear(VALUE self) {
  VALUE result = Qtrue;

  kclear_c();
  check_spice_error();

  return result;
}
/*!@fn VALUE ktotal(int argc, VALUE *argv, VALUE self)
 * @brief  Return the current number of kernels that have been loaded 
 *    via the KEEPER interface that are of a specified type. 
 * @brief http://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/ktotal_c.html
 * @param NONE
 *
 * @return a fixnum representing the number of kernels
 */
VALUE ktotal(int argc, VALUE *argv, VALUE self) {
  VALUE result = Qnil;
  SpiceInt numkernels;

  if (argc > 0) {
    rb_raise(rb_eArgError, "no parameters here");
    return result;
  }

  ktotal_c("ALL", &numkernels);
  check_spice_error();

  return INT2FIX(numkernels);
}
/*!@fn VALUE str2et(int argc, VALUE *argv, VALUE self)
 * @brief Convert a string representing an epoch to a double precision
 *    value representing the number of TDB seconds past the J2000
 *    epoch corresponding to the input epoch.
 * @brief http://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/str2et_c.html
 *
 *  @param a time string 
 *  @returns a float representing et
 *
 */
VALUE str2et(int argc, VALUE *argv, VALUE self) {
  double et;

  if (argc != 1) {
    rb_raise(rb_eArgError, "need one parameter, the time string");
    return Qnil;
  }

  Check_Type(argv[0], T_STRING);

  str2et_c(StringValuePtr(argv[0]), &et);

  check_spice_error();

  return rb_float_new(et);
}
/*!@fn VALUE subpnt(int argc, VALUE *argv, VALUE self)
 *@brief  Compute the rectangular coordinates of the sub-observer point on 
 *   a target body at a specified epoch, optionally corrected for 
 *     light time and stellar aberration. 
 *
 *@brief http://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/subpnt_c.html
 *
 *  @param method    
 *     Computation method. EX: "Near point: ellipsoid" or  "Intercept: ellipsoid"
 *  @param target    
 *     Name of target body. EX: "MOON"
 *  @param et     
 *     Epoch in ephemeris seconds past J2000 TDB. 
 *  @param fixref  
 *      Body-fixed, body-centered target body frame. EX: "iau_moon"
 *  @param abcorr   
 *       Aberration correction. Please see @ref abcorr
 *  @param obsrvr   
 *       Name of observing body. EX: "LRO" 
 *
 *@returns  a ruby array [[spointX, spointY, spointZ], trgepc, [srfvecX, srfvecY, srfvecZ]]
 *@returns  spoint (vector) = Sub-observer point on the target body. 
 *@returns  trgepc          = Sub-observer point epoch. 
 *@returns  srfvec (vector) = Vector from observer to sub-observer point. 
 *
 */
VALUE subpnt(int argc, VALUE *argv, VALUE self) {
  VALUE result = Qnil;
  double spoint[3];
  double trgepc;
  double srfvec[3];
  VALUE rb_spoint, rb_srfvec;
  SpiceInt i;

  if (argc != 6) {
    rb_raise(rb_eArgError, "need 6 parameters!");
    return Qnil;
  }

  Check_Type(argv[0], T_STRING);
  Check_Type(argv[1], T_STRING);
  Check_Type(argv[2], T_FLOAT);
  Check_Type(argv[3], T_STRING);
  Check_Type(argv[4], T_STRING);
  Check_Type(argv[5], T_STRING);

  subpnt_c(StringValuePtr(argv[0]), StringValuePtr(argv[1]), NUM2DBL(argv[2]), StringValuePtr(argv[3]), StringValuePtr(argv[4]), StringValuePtr(argv[5]), spoint, &trgepc, srfvec);

  rb_spoint = rb_ary_new();
  for (i=0; i < 3; i++)
    rb_ary_push(rb_spoint, rb_float_new(spoint[i]));

  rb_srfvec = rb_ary_new();
  for (i=0; i < 3; i++)
    rb_ary_push(rb_srfvec, rb_float_new(srfvec[i]));

  result = rb_ary_new();
  rb_ary_push(result, rb_spoint);
  rb_ary_push(result, rb_float_new(trgepc));
  rb_ary_push(result, rb_srfvec);

  check_spice_error();

  return result;
}
/*!@fn VALUE reclat(int argc, VALUE* argv, VALUE self)
* @brief Convert from rectangular coordinates to latitudinal coordinates.
* @brief http://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/reclat_c.html
*
*@param rectan
*       The rectangular coordinates of the input point.  `rectan'
*       is a 3-vector.
*
*@returns [radius, lon, lat]
* @returns returns a ruby array that can be saved as follows
* @code radius, lon, lat = reclat(rectan) @endcode
*/
VALUE reclat(int argc, VALUE* argv, VALUE self) {
  double spoint[3];
  SpiceInt i;
  double radius, lon, lat;

  if (argc != 1) {
    rb_raise(rb_eArgError, "Need just 1 parameter");
    return Qnil;
  }

  Check_Type(argv[0], T_ARRAY);

  if (RARRAY_LEN(argv[0]) != 3) {
    rb_raise(rb_eArgError, "The array should have 3 items in it");
    return Qnil;
  }

  for (i=0; i < 3; i++) {
    Check_Type(RARRAY_PTR(argv[0])[i], T_FLOAT);
    spoint[i] = NUM2DBL(RARRAY_PTR(argv[0])[i]);
  }

  reclat_c(spoint, &radius, &lon, &lat);

  check_spice_error();

  return rb_ary_new3(3, rb_float_new(radius), rb_float_new(lon), rb_float_new(lat));
}
/*!@fn VALUE dpr(int argc, VALUE *argv, VALUE self)
*@brief Return the number of degrees per radian.
*@brief http://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/dpr_c.html
*
*@param NONE
*
*@returns 180.0/acos(-1)
*
*EX: 
* @code angleDEG = angleRAD*dpr() @endcode
*
*/
VALUE dpr(int argc, VALUE *argv, VALUE self) {
  if (argc != 0) {
    rb_raise(rb_eArgError, "no args, go away!");
    return Qnil;
  }

  return rb_float_new(dpr_c());
}

/*!@fn VALUE eul2m(int argc, VALUE *argv, VALUE self)
*@brief Construct a rotation matrix from a set of Euler angles.
*@brief http://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/eul2m_c.html
*
*@param  angle3 (IN RADIANS!)
*@param  angle2 (IN RADIANS!)
*@param  angle1 (IN RADIANS!)
*@param  axis3  Axis number of the third rotation axes.
*@param  axis2  Axis number of the second rotation axes.
*@param  axis1  Axis number of the first rotation axes.
*
*@returns A 3X3 ruby array representing the product of the 3 rotations.
*
*/
VALUE eul2m(int argc, VALUE *argv, VALUE self){
	double angle1, angle2, angle3;
	int axis1, axis2, axis3;
	int i;
	double r[3][3];

	if (argc != 6){
		rb_raise(rb_eArgError, "Need 6 Parameters!");
		return Qnil;
	}
	for( i = 0; i < 3; i++)
		Check_Type(argv[i], T_FLOAT);
	for( i = 3; i < 6; i++)
		Check_Type(argv[i], T_FIXNUM);

	angle3 = NUM2DBL(argv[0]);
	angle2 = NUM2DBL(argv[1]);
	angle1 = NUM2DBL(argv[2]);
	axis3  = FIX2INT(argv[3]);
	axis2  = FIX2INT(argv[4]);
	axis1  = FIX2INT(argv[5]);

	eul2m_c(angle3, angle2,  angle1, axis3, axis2, axis1, r);

	check_spice_error();

	return rb_ary_new3(3, rb_ary_new3(3, rb_float_new(r[0][0]), rb_float_new(r[0][1]), rb_float_new(r[0][2])),  rb_ary_new3(3, rb_float_new(r[1][0]), rb_float_new(r[1][1]), rb_float_new(r[1][2])), rb_ary_new3(3, rb_float_new(r[2][0]), rb_float_new(r[2][1]), rb_float_new(r[2][2])));

}

/*!@fn VALUE vnorm(int argc, VALUE *argv, VALUE self)
*@brief  Compute the magnitude of a double precision, 3-dimensional vector.
*@brief http://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/vnorm_c.html
*
*@param v1  
*        Vector whose magnitude is to be found.
*
*@returns A ruby float representing the magnitude of the vector
*
*EX:
* @code vector = 1,2,3
* mag = vnorm(vector)@endcode
*/
VALUE vnorm(int argc, VALUE *argv, VALUE self) {
  double srfvec[3], res;
  SpiceInt i;

  if (argc != 1) {
    rb_raise(rb_eArgError, "Need just 1 parameter");
    return Qnil;
  }

  Check_Type(argv[0], T_ARRAY);

  if (RARRAY_LEN(argv[0]) != 3) {
    rb_raise(rb_eArgError, "The array should have 3 items in it");
    return Qnil;
  }

  for (i=0; i < 3; i++) {
    Check_Type(RARRAY_PTR(argv[0])[i], T_FLOAT);
    srfvec[i] = NUM2DBL(RARRAY_PTR(argv[0])[i]);
  }

  res = vnorm_c(srfvec);

  check_spice_error();

  return rb_float_new(res);
}

/*!@fn VALUE vperp(VALUE self, VALUE v1, VALUE v2)
*@brief  Find the component of a vector that is perpendicular to a second vector. 
*@brief http://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/vperp_c.html
*
*@param v1  
*        The vector whose orthogonal component is sought.
*@param v2  
*        The vector used as the orthogonal reference.
*
*@returns A ruby aray[3] representing The component of a orthogonal to b.
*
*/
VALUE vperp(VALUE self, VALUE v1, VALUE v2) {
  double vec1[3], vec2[3], result[3];
  SpiceInt i;

  Check_Type(v1, T_ARRAY);
  Check_Type(v2, T_ARRAY);

  if (RARRAY_LEN(v1) != 3) {
    rb_raise(rb_eArgError, "The array should have 3 items in it");
    return Qnil;
  }
  if (RARRAY_LEN(v2) != 3) {
    rb_raise(rb_eArgError, "The array should have 3 items in it");
    return Qnil;
  }

  for (i=0; i < 3; i++) {
    Check_Type(RARRAY_PTR(v1)[i], T_FLOAT);
    vec1[i] = NUM2DBL(RARRAY_PTR(v1)[i]);
  }
  for (i=0; i < 3; i++) {
    Check_Type(RARRAY_PTR(v2)[i], T_FLOAT);
    vec2[i] = NUM2DBL(RARRAY_PTR(v2)[i]);
  }

  vperp_c(vec1, vec2, result);

  check_spice_error();

  return rb_ary_new3(3, rb_float_new(result[0]), rb_float_new(result[1]), rb_float_new(result[2]));
}

/*!@fn VALUE ilumin(int argc, VALUE *argv, VALUE self)
* @brief Find the illumination angles (phase, solar incidence, and
*   emission) at a specified surface point of a target body.
* @brief http://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/ilumin_c.html
*
*@param method 
*Computation method. 
*The only choice currently supported is "Ellipsoid"
*@param target  
*Name of target body. EX: "MOON" 
*@param et     
*Epoch in ephemeris seconds past J2000 TDB.
*@param fixref   
*Body-fixed, body-centered target body frame. EX "iau_moon"
*@param abcorr     
*Desired aberration correction. Please see @ref abcorr
*@param obsrvr   
*Name of observing body. EX: "LRO"
*@param spoint    
*Body-fixed coordinates of a target surface point. a vector
*
*
*@returns a ruby array:
*@returns [trgepc, [srfvecX, srfvecY, srfvecZ], phaseAngle, solarAngle, emissionAngle]
*@returns NOTE: all angles are in radians
*
*/
VALUE ilumin(int argc, VALUE *argv, VALUE self) {
  double spoint[3];
  double trgepc;
  double srfvec[3];
  double phase, solar, emissn;
  VALUE rb_srfvec;
  SpiceInt i;

  if (argc != 7) {
    rb_raise(rb_eArgError, "need 7 parameters!");
    return Qnil;
  }

  Check_Type(argv[0], T_STRING);
  Check_Type(argv[1], T_STRING);
  Check_Type(argv[2], T_FLOAT);
  Check_Type(argv[3], T_STRING);
  Check_Type(argv[4], T_STRING);
  Check_Type(argv[5], T_STRING);
  Check_Type(argv[6], T_ARRAY);

  if (RARRAY_LEN(argv[6]) != 3) {
    rb_raise(rb_eArgError, "The array should have 3 items in it");
    return Qnil;
  }

  for (i=0; i < 3; i++) {
    Check_Type(RARRAY_PTR(argv[6])[i], T_FLOAT);
    spoint[i] = NUM2DBL(RARRAY_PTR(argv[6])[i]);
  }

  ilumin_c(StringValuePtr(argv[0]), StringValuePtr(argv[1]), NUM2DBL(argv[2]), StringValuePtr(argv[3]), StringValuePtr(argv[4]), StringValuePtr(argv[5]), spoint, &trgepc, srfvec, &phase, &solar, &emissn);

  rb_srfvec = rb_ary_new();
  for (i=0; i < 3; i++)
    rb_ary_push(rb_srfvec, rb_float_new(srfvec[i]));

  check_spice_error();

  return rb_ary_new3(5, rb_float_new(trgepc), rb_srfvec, rb_float_new(phase), rb_float_new(solar), rb_float_new(emissn));
}
/*!@fn VALUE spkcov(int argc, VALUE *argv, VALUE self)
*@brief Find the coverage window for a specified ephemeris object in a 
*   specified SPK file.
*@brief http://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/spkcov_c.html
*
*@param spk 
*       name of the spk file, Can be a single spk OR an array of spks
*@param idcode
*	is the integer ID code of an object for which 
*	ephemeris data are expected to exist in the 
*	specified SPK file.
*
*@returns An array of start/end time pairs which the spk covers 
*@returns EX: [[start,end], [start_end], ...]
*
*
*/
VALUE spkcov(int argc, VALUE *argv, VALUE self) {
  SPICEDOUBLE_CELL (cover, 2000);
  SpiceInt i;
  SpiceInt niv;
  double b, e;
  VALUE result;

  if (argc != 2) {
    rb_raise(rb_eArgError, "We need two args!");
    return Qnil;
  }

 
  if (TYPE(argv[0]) != T_STRING && TYPE(argv[0]) != T_ARRAY) {
    rb_raise(rb_eArgError, "gimme some spks!");
    return Qnil;
  }

  Check_Type(argv[1], T_FIXNUM);

  scard_c(0, &cover);

  if (TYPE(argv[0]) == T_STRING) {
    spkcov_c(StringValuePtr(argv[0]), FIX2INT(argv[1]), &cover);
  } else {
    for (i=0; i < RARRAY_LEN(argv[0]); i++)
      spkcov_c(StringValuePtr(RARRAY_PTR(argv[0])[i]), FIX2INT(argv[1]), &cover);
  }

  niv = wncard_c(&cover);

  if (niv == 0)
    return Qnil;

  result = rb_ary_new2(niv);

  for (i = 0; i < niv; i++) {
    wnfetd_c(&cover, i, &b, &e);
    printf("%f %f\n",b, e);
    rb_ary_push(result, rb_ary_new3(2, INT2FIX(b), INT2FIX(e)));
  }
  /*build up an array of arrays with start/end time pairs*/

  check_spice_error();

  return result;
}
/*!@fn VALUE bodn2c(int argc, VALUE *argv, VALUE self)
*@brief Translate the name of a body or object to the corresponding SPICE integer ID code.
*@brief http://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/bodn2c_c.html
*
*@param name 
	the name of a body or object, such as a planet,
        satellite, comet, asteroid, barycenter, DSN station,
        spacecraft, or instrument, that is "known" to the SPICE
        system, whether through hard-coded registration or
        run-time registration in the SPICE kernel pool.
        Case and leading and trailing blanks in `name'
        are not significant.  However when a name is made
        up of more than one word, they must be separated by
        at least one blank. \n
	  That is, all of the following strings are equivalent names: \n 
                "JUPITER BARYCENTER" \n
                "Jupiter Barycenter" \n
                "JUPITER BARYCENTER   " \n
                "JUPITER    BARYCENTER" \n 
                "   JUPITER BARYCENTER" \n 
         However, "JUPITERBARYCENTER" is not equivalent to
         the names above.
*
*
*@returns nil           IF no code is found
*@returns code(float)  IF a code is found
*/
VALUE bodn2c(int argc, VALUE *argv, VALUE self) {
  SpiceBoolean found;
  SpiceInt code;

  if (argc != 1) {
    rb_raise(rb_eArgError, "need just one arg please!");
    return Qnil;
  }

  Check_Type(argv[0], T_STRING);

  bodn2c_c(StringValuePtr(argv[0]), &code, &found);

  if (found == SPICEFALSE)
    return Qnil;
  
  check_spice_error();

  return INT2FIX(code);
}
/*!@fn VALUE et2rb(VALUE self, VALUE arg)
*@brief converts et into a ruby DateTime object
*
*@param et 
*	a time value in et
*
*@returns A ruby DateTime object
*/
VALUE et2rb(VALUE self, VALUE arg) {
  double et;
  char time_str[64];
  char tmpstr[1000];
  VALUE datetime;

  et = NUM2DBL(arg);

  et2utc_c(et, "ISOC", 3, 64, time_str);

  rb_eval_string("require 'date'");
  sprintf(tmpstr, "dt = DateTime.parse(\"%s\"); Time.utc(dt.year, dt.month, dt.day, dt.hour, dt.min, dt.sec, (dt.sec_fraction*(RUBY_VERSION =~ /^1\\.8\\./ ? 24*60*60 : 1)*1000000).to_i)", time_str);
  datetime = rb_eval_string(tmpstr);

  check_spice_error();

  return datetime;
}
/*!@fn VALUE rb2et(VALUE self, VALUE arg)
*@brief convert a ruby DateTime object into et
*
*@param DateTime 
*	A ruby Date Time Object
*
*@returns et
*/
VALUE rb2et(VALUE self, VALUE arg) {
  double et;
  VALUE rb_tmp1, rb_tmp2;
  int usec;
  char tmpstr[1000];

  usec = FIX2INT(rb_funcall(arg, rb_intern("usec"), 0));
  rb_tmp1 = rb_funcall(arg, rb_intern("utc"), 0);

  rb_tmp2 = rb_funcall(rb_tmp1, rb_intern("strftime"), 1, rb_str_new2("%Y-%m-%dT%H:%M:%S"));

  sprintf(tmpstr, "%s.%06d", StringValuePtr(rb_tmp2), usec);

  str2et_c(tmpstr, &et);

  check_spice_error();

  return rb_float_new(et);
}
/*!@fn VALUE subslr(int argc, VALUE *argv, VALUE self)
*@brief Compute the rectangular coordinates of the sub-solar point on 
* 	 a target body at a specified epoch, optionally corrected for 
*	 light time and stellar aberration. 
*@brief http://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/subslr_c.html
*
*@param method    
*   Computation method. EX: "Near point: ellipsoid" or  "Intercept: ellipsoid"
*@param target    
*   Name of target body. EX: "MOON"
*@param et     
*   Epoch in ephemeris seconds past J2000 TDB. 
*@param fixref  
*    Body-fixed, body-centered target body frame. EX: "iau_moon"
*@param abcorr   
*     Aberration correction. Please see @ref abcorr
*@param obsrvr   
*     Name of observing body. EX: "LRO" 
*
*@returns  a ruby array [[spointX, spointY, spointZ], trgepc, [srfvecX, srfvecY, srfvecZ]]
*@returns  spoint (vector) = Sub-observer point on the target body. 
*@returns  trgepc          = Sub-observer point epoch. 
*@returns  srfvec (vector) = Vector from observer to sub-solar point. 
*
*/
VALUE subslr(int argc, VALUE *argv, VALUE self) {
  VALUE result = Qnil;
  double spoint[3];
  double trgepc;
  double srfvec[3];
  VALUE rb_spoint, rb_srfvec;
  SpiceInt i;

  if (argc != 6) {
    rb_raise(rb_eArgError, "need 6 parameters!");
    return Qnil;
  }

  Check_Type(argv[0], T_STRING);
  Check_Type(argv[1], T_STRING);
  Check_Type(argv[2], T_FLOAT);
  Check_Type(argv[3], T_STRING);
  Check_Type(argv[4], T_STRING);
  Check_Type(argv[5], T_STRING);

  subslr_c(StringValuePtr(argv[0]), StringValuePtr(argv[1]), NUM2DBL(argv[2]), StringValuePtr(argv[3]), StringValuePtr(argv[4]), StringValuePtr(argv[5]), spoint, &trgepc, srfvec);

  rb_spoint = rb_ary_new();
  for (i=0; i < 3; i++)
    rb_ary_push(rb_spoint, rb_float_new(spoint[i]));

  rb_srfvec = rb_ary_new();
  for (i=0; i < 3; i++)
    rb_ary_push(rb_srfvec, rb_float_new(srfvec[i]));

  result = rb_ary_new();
  rb_ary_push(result, rb_spoint);
  rb_ary_push(result, rb_float_new(trgepc));
  rb_ary_push(result, rb_srfvec);

  check_spice_error();

  return result;
}
/*!@fn VALUE sincpt(int argc, VALUE *argv, VALUE self)
*@brief Given an observer and a direction vector defining a ray, compute 
*   the surface intercept of the ray on a target body at a specified 
*   epoch, optionally corrected for light time and stellar 
*   aberration. 
*@brief http://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/sincpt_c.html
*
*@param   method        Computation method. 
*@param   target        Name of target body. 
*@param   et            Epoch in ephemeris seconds past J2000 TDB. 
*@param   fixref        Body-fixed, body-centered target body frame. 
*@param   abcorr        Aberration correction. please see @ref abcorr
*@param   obsrvr        Name of observing body. 
*@param   dref          Reference frame of ray's direction vector. 
*@param   dvec          Ray's direction vector. 
*
*@returns nil if no intersection is found. Otherwise:
*@returns A ruby array
*@returns [[spointX, spointY, spointZ],trgepc,[srfvecX, srfvecY, srfvecZ]]
*
*/
VALUE sincpt(int argc, VALUE *argv, VALUE self) {
  double spoint[3] = {0.0,0.0,0.0};
  double trgepc;
  double srfvec[3];
  double dvec[3];
  SpiceBoolean found;
  VALUE rb_srfvec, rb_spoint;
  SpiceInt i;

  if (argc != 8) {
    rb_raise(rb_eArgError, "need 7 parameters!");
    return Qnil;
  }

  Check_Type(argv[0], T_STRING);
  Check_Type(argv[1], T_STRING);
  Check_Type(argv[2], T_FLOAT);
  Check_Type(argv[3], T_STRING);
  Check_Type(argv[4], T_STRING);
  Check_Type(argv[5], T_STRING);
  Check_Type(argv[6], T_STRING);
  Check_Type(argv[7], T_ARRAY);

  if (RARRAY_LEN(argv[7]) != 3) {
    rb_raise(rb_eArgError, "The array should have 3 items in it");
    return Qnil;
  }

  for (i=0; i < 3; i++) {
    Check_Type(RARRAY_PTR(argv[7])[i], T_FLOAT);
    dvec[i] = NUM2DBL(RARRAY_PTR(argv[7])[i]);
  }

  sincpt_c(StringValuePtr(argv[0]), StringValuePtr(argv[1]), NUM2DBL(argv[2]), StringValuePtr(argv[3]), StringValuePtr(argv[4]), StringValuePtr(argv[5]), StringValuePtr(argv[6]), dvec, spoint, &trgepc, srfvec, &found);

  if (found == 0) {
    rb_raise(rb_eRuntimeError, "Intercept point not found");
    return Qnil;
  }

  rb_spoint = rb_ary_new();
  for (i=0; i < 3; i++)
    rb_ary_push(rb_spoint, rb_float_new(spoint[i]));

  rb_srfvec = rb_ary_new();
  for (i=0; i < 3; i++)
    rb_ary_push(rb_srfvec, rb_float_new(srfvec[i]));

  check_spice_error();

  return rb_ary_new3(3, rb_spoint, rb_float_new(trgepc), rb_srfvec);
}
/*!@fn VALUE getfov(VALUE self, VALUE arg)
*@brief  This routine returns the field-of-view (FOV) parameters for a
*   specified instrument.
*@brief http://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/getfov_c.html
*
*@param   instid        NAIF ID of an instrument.
*
*
*@returns A ruby Array:
*@returns [shape, frame, [bsight], [bounds]]
*@returns   shape:         Instrument FOV shape. 
*@returns   frame:         Name of the frame in which FOV vectors are defined. 
*@returns   bsight:        Boresight vector. 
*@returns   bounds:        FOV boundary vectors. 
*
*/
VALUE getfov(VALUE self, VALUE arg) {
  char shape[32];
  char frame[32];
  double bsight[3];
  SpiceInt n;
  double bounds[10][3];

  SpiceInt i;
  VALUE rb_bsight, rb_bounds;

  Check_Type(arg, T_FIXNUM);

  getfov_c(FIX2INT(arg), 10, 32, 32, shape, frame, bsight, &n, bounds);

  rb_bsight = rb_ary_new();
  for (i=0; i < 3; i++)
    rb_ary_push(rb_bsight, rb_float_new(bsight[i]));

  rb_bounds = rb_ary_new2(n);
  for (i = 0; i < n; i++)
    rb_ary_push(rb_bounds, rb_ary_new3(3, rb_float_new(bounds[i][0]), rb_float_new(bounds[i][1]), rb_float_new(bounds[i][2])));

  check_spice_error();

  return rb_ary_new3(4, rb_str_new2(shape), rb_str_new2(frame), rb_bsight, rb_bounds);
}
/*!@fn VALUE gfdist(int argc, VALUE *argv, VALUE self)
*@brief  Return the time window over which a specified constraint on   
* observer-target distance is met. 
*@brief http://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/gfdist_c.html
*
*@param   target               Name of the target body. 
*@param   abcorr               Aberration correction flag. 
*@param   obsrvr               Name of the observing body. 
*@param   relate               Relational operator.\n
Acceptable value's are:\n
	 ">"      Distance is greater than the reference value `refval'. \n
	"="      Distance is equal to the reference value `refval'. \n
	"<"      Distance is less than the reference value `refval'. \n
	"ABSMAX"  Distance is at an absolute maximum. \n
	"ABSMIN"  Distance is at an absolute  minimum. \n
	"LOCMAX"  Distance is at a local maximum. \n
	"LOCMIN"  Distance is at a local minimum. \n\n
*@param   refval               Reference value is the reference value used together with the argument `relate' to define an equality
*@param   adjust               Adjustment value for absolute extrema searches. 
*	        adjust a parameter used to modify searches for absolute
*               extrema: when `relate' is set to "ABSMAX" or "ABSMIN"
*               and `adjust' is set to a positive value, gfdist_c will
*               find times when the observer-target distance is within
*               `adjust' km of the specified extreme value.\n
*               If `adjust' is non-zero and a search for an absolute
*               minimum `min' is performed, the result window contains
*               time intervals when the observer-target distance has
*               values between `min' and min+adjust.\n
*               If the search is for an absolute maximum `max', the
*               corresponding range is from max-adjust to `max'.\n
*               `adjust' is not used for searches for local extrema,
*               equality or inequality conditions
*@param   step                 Step size used for locating extrema and roots. 
*@param   nintvls              Workspace window interval count. (number of intervals to return 
*@param   cnfine **            SPICE window to which the search is confined. 
* cfntime is a spiceDouble window. It must be initalized with wninsd()
*@param etstart   Instead of a cfn window you can specify an et start and stop time (OSX)
*@param etstop
*
*@returns A 2-demesional ruby array of windows
*@returns [[start,end],[start,end]...]
*@returns times are returned in et
*/
VALUE gfdist(int argc, VALUE *argv, VALUE self) {
	SpiceInt i, niv;
	double beg, end;
	VALUE result;
	SPICEDOUBLE_CELL(intv, 5000);
	SPICEDOUBLE_CELL(win, 5000);
	SpiceCell *ptr;

	if (argc < 9 || argc > 10){
		rb_raise(rb_eArgError, "need 9 parameters!");
		return Qnil;
	}

	Check_Type(argv[0], T_STRING);/*target*/
	Check_Type(argv[1], T_STRING);/*abcorr*/
	Check_Type(argv[2], T_STRING);/*obsrvr*/
	Check_Type(argv[3], T_STRING);/*relate*/
	Check_Type(argv[4], T_FLOAT);/*refval*/
	Check_Type(argv[5], T_FLOAT);/*adj*/
	Check_Type(argv[6], T_FLOAT);/*step*/
	Check_Type(argv[7], T_FIXNUM);/*nintvls*/


	if (argc == 9){
		Check_Type(argv[8], T_DATA);/*spicedouble_Cell window*/

		Data_Get_Struct(argv[8], SpiceCell, ptr);
		gfdist_c (StringValuePtr(argv[0]), 
				StringValuePtr(argv[1]), 
				StringValuePtr(argv[2]), 
				StringValuePtr(argv[3]), 
				NUM2DBL(argv[4]), 
				NUM2DBL(argv[5]), 
				NUM2DBL(argv[6]), 
				FIX2INT(argv[7]),
				ptr , 
				&intv );
	}
	else{
		Check_Type(argv[8], T_FLOAT);/*start et*/
		Check_Type(argv[9], T_FLOAT);/*stop et*/

		wninsd_c(NUM2DBL(argv[8]), NUM2DBL(argv[9]), &win);
		gfdist_c(StringValuePtr(argv[0]), 
				StringValuePtr(argv[1]), 
				StringValuePtr(argv[2]), 
				StringValuePtr(argv[3]), 
				NUM2DBL(argv[4]), 
				NUM2DBL(argv[5]), 
				NUM2DBL(argv[6]), 
				FIX2INT(argv[7]),
				&win , 
				&intv );
	}
	check_spice_error();

	niv = wncard_c(&intv);
	if (niv == 0)
		return Qnil;

	result = rb_ary_new2(niv);

	for (i = 0; i < niv; i++) {
		wnfetd_c(&intv, i, &beg, &end);
		rb_ary_push(result, rb_ary_new3(2, INT2FIX(beg), INT2FIX(end)));
	}

	return result;
}

/*!@fn VALUE gfsntc(int argc, VALUE *argv, VALUE self)
*@brief    Determine time intervals for which a coordinate of an
*          surface intercept position vector satisfies a numerical constraint.
*
*@brief http://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/gfsntc_c.html
*
*@param   targ 	  Name of the target body  
*@param   fixref  Body fixed frame associated with 'target'
*@param   method  Name of method type for surface intercept calculation
*@param   abcorr  Aberration correction flag  
*@param   obsrvr  Name of the observing body
*@param   dref    Reference frame of direction vector 'dvec'
*@param   dvec    Pointing direction vector from 'obsrvr' (a ruby 3 vector of floats)
*@param   crdsys  Name of the coordinate system containing COORD
*@param   coord   Name of the coordinate of interest
*@param   relate  Operator that either looks for an extreme value
		  (max, min, local, absolute) or compares the
		  coordinate value and refval
*@param   refval  Reference value
*@param   adjust  Adjustment value for absolute extrema searches
*@param   step    Step size used for locating extrema and roots
*@param   nintvls Workspace window interval count
*@param   cnfine        SPICE window to which the search is restricted
* cfntime is a spiceDouble window. It must be initalized with wninsd() 
* or a start_et and end_et can be used
*
*
*@returns A 2-demesional ruby array of windows
*@returns [[start,end],[start,end]...]
*@returns times are returned in et
*
*/
VALUE gfsntc(int argc, VALUE *argv, VALUE self) {
	double dvec[3];
	SpiceInt i, niv;
	double beg, end;
	VALUE result;
	SPICEDOUBLE_CELL(intv, 5000);
	SPICEDOUBLE_CELL(win, 5000);
	SpiceCell *ptr;
	if (argc < 15  || argc > 16){
		rb_raise(rb_eArgError, "need 15 parameters!");
		return Qnil;
	}

	Check_Type(argv[0], T_STRING);/*//target*/
	Check_Type(argv[1], T_STRING);/*fixref*/
	Check_Type(argv[2], T_STRING);/*method*/
	Check_Type(argv[3], T_STRING);/*abcorr*/
	Check_Type(argv[4], T_STRING);/*obsrvr*/
	Check_Type(argv[5], T_STRING);/*dref*/
	Check_Type(argv[6], T_ARRAY); /*dvec*/
	Check_Type(argv[7], T_STRING);/*crdsys*/
	Check_Type(argv[8], T_STRING);/*coord*/
	Check_Type(argv[9], T_STRING);/*relate*/
	Check_Type(argv[10], T_FLOAT); /*refval*/
	Check_Type(argv[11], T_FLOAT); /*adjust*/
	Check_Type(argv[12], T_FLOAT); /*step*/
	Check_Type(argv[13], T_FIXNUM);/*nintvls*/

	for (i=0; i < 3; i++){
		Check_Type(RARRAY_PTR(argv[6])[i], T_FLOAT);
		dvec[i] = NUM2DBL(RARRAY_PTR(argv[6])[i]);
	}

	if (argc == 15){
		Check_Type(argv[14], T_DATA);/*spicedouble_Cell window*/
		Data_Get_Struct(argv[14], SpiceCell, ptr);

		gfsntc_c(StringValuePtr(argv[0]),
				StringValuePtr(argv[1]),
				StringValuePtr(argv[2]),
				StringValuePtr(argv[3]),
				StringValuePtr(argv[4]),
				StringValuePtr(argv[5]),
				dvec,
				StringValuePtr(argv[7]),
				StringValuePtr(argv[8]),
				StringValuePtr(argv[9]),
				NUM2DBL(argv[10]),
				NUM2DBL(argv[11]),
				NUM2DBL(argv[12]),
				FIX2INT(argv[13]),
				ptr,
				&intv);
	}
	else {
		Check_Type(argv[14], T_FLOAT);/*start et*/
		Check_Type(argv[15], T_FLOAT);/*stop et*/
		wninsd_c(NUM2DBL(argv[14]), NUM2DBL(argv[15]), &win);

		gfsntc_c ( StringValuePtr(argv[0]),
				StringValuePtr(argv[1]),
				StringValuePtr(argv[2]),
				StringValuePtr(argv[3]),
				StringValuePtr(argv[4]),
				StringValuePtr(argv[5]),
				dvec,
				StringValuePtr(argv[7]),
				StringValuePtr(argv[8]),
				StringValuePtr(argv[9]),
				NUM2DBL(argv[10]),
				NUM2DBL(argv[11]),
				NUM2DBL(argv[12]),
				FIX2INT(argv[13]),
				&win,
				&intv);

	}
	check_spice_error();

	niv = wncard_c(&intv);
	if (niv == 0)
		return Qnil;

	result = rb_ary_new2(niv);

	for (i = 0; i < niv; i++) {
		wnfetd_c(&intv, i, &beg, &end);
		rb_ary_push(result, rb_ary_new3(2, rb_float_new(beg), rb_float_new(end)));
	}

	return result;
}

/*!@fn VALUE gfsep(int argc, VALUE *argv, VALUE self)
*@brief Determine time intervals when the angular separation between
*   the position vectors of two target bodies relative to an observer
*   satisfies a numerical relationship.
*@brief http://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/gfsep_c.html
*
*@param   targ1         Name of first body
*@param   shape1        Name of shape model describing the first body
*@param   frame1        The body-fixed reference frame of the first body
*@param   targ2         Name of second body
*@param   shape2        Name of the shape model describing the second body
*@param   frame2        The body-fixed reference frame of the second body
*@param   abcorr        Aberration correction flag. See @ref abcorr
*@param   obsrvr        Name of the observing body
*@param   relate        Operator that either looks for an extreme value
*                  (max, min, local, absolute) or compares the
*                  angular separation value and refval
*@param   refval        Reference value
*@param   adjust        Absolute extremum adjustment value
*@param   step          Step size in seconds for finding angular separation events
*@param   nintvls       Workspace window interval count
*@param   cnfine        SPICE window to which the search is restricted
* cfntime is a spiceDouble window. It must be initalized with wninsd()
*
*@returns A 2-demesional ruby array of windows
*@returns [[start,end],[start,end]...]
*@returns times are returned in et
*
*/

VALUE gfsep(int argc, VALUE *argv, VALUE self) {
	SpiceInt i, niv;
	double beg, end;
	VALUE result;
	SPICEDOUBLE_CELL(intv, 5000);
	SpiceCell *ptr;
	SPICEDOUBLE_CELL(win, 5000);

	if (argc < 14 || argc >> 15){
		rb_raise(rb_eArgError, "need 14 parameters!");
		return Qnil;
	}
	Check_Type(argv[0], T_STRING);/*target1*/
	Check_Type(argv[1], T_STRING);/*shape1*/
	Check_Type(argv[2], T_STRING);/*frame1*/
	Check_Type(argv[3], T_STRING);/*targ2*/
	Check_Type(argv[4], T_STRING);/*shape2*/
	Check_Type(argv[5], T_STRING);/*frame2*/
	Check_Type(argv[6], T_STRING);/*abcorr*/
	Check_Type(argv[7], T_STRING);/*obsvr*/
	Check_Type(argv[8], T_STRING);/*relate*/
	Check_Type(argv[9], T_FLOAT);/*refval*/
	Check_Type(argv[10], T_FLOAT);/*adjust*/
	Check_Type(argv[11], T_FLOAT);/*step*/
	Check_Type(argv[12], T_FIXNUM);/*nintvls*/


	if (argc == 14){
		Check_Type(argv[13], T_DATA);/*cfnine*/
		Data_Get_Struct(argv[13], SpiceCell, ptr);

		gfsep_c(StringValuePtr(argv[0]),
				StringValuePtr(argv[1]),
				StringValuePtr(argv[2]),
				StringValuePtr(argv[3]),
				StringValuePtr(argv[4]),
				StringValuePtr(argv[5]),
				StringValuePtr(argv[6]),
				StringValuePtr(argv[7]),
				StringValuePtr(argv[8]),
				NUM2DBL(argv[9]),
				NUM2DBL(argv[10]),
				NUM2DBL(argv[11]),
				FIX2INT(argv[12]),
				ptr,
				&intv);
	}
	else{
		Check_Type(argv[13], T_FLOAT);/*start et*/
		Check_Type(argv[14], T_FLOAT);/*stop et*/

		wninsd_c(NUM2DBL(argv[13]), NUM2DBL(argv[14]), &win);

		gfsep_c(StringValuePtr(argv[0]),
				StringValuePtr(argv[1]),
				StringValuePtr(argv[2]),
				StringValuePtr(argv[3]),
				StringValuePtr(argv[4]),
				StringValuePtr(argv[5]),
				StringValuePtr(argv[6]),
				StringValuePtr(argv[7]),
				StringValuePtr(argv[8]),
				NUM2DBL(argv[9]),
				NUM2DBL(argv[10]),
				NUM2DBL(argv[11]),
				FIX2INT(argv[12]),
				&win,
				&intv);

	}


	check_spice_error();
	niv = wncard_c(&intv);
	if (niv == 0)
		return Qnil;

	result = rb_ary_new2(niv);

	for (i = 0; i < niv; i++){
		wnfetd_c(&intv, i, &beg, &end);
		rb_ary_push(result, rb_ary_new3(2, INT2FIX(beg), INT2FIX(end)));
	}

	return result;

}
/*!@fn VALUE gftfov(int argc, VALUE *argv, VALUE self)
*@brief Determine time intervals when a specified ephemeris object 
*   intersects the space bounded by the field-of-view (FOV) of a 
*   specified instrument. 
*@brief http://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/gftfov_c.html
*
*@param   inst                 Name of the instrument. 
*@param   target               Name of the target body. 
*@param   tshape               Type of shape model used for target body. 
*@param   tframe               Body-fixed, body-centered frame for target body. 
*@param   abcorr               Aberration correction flag. see @ref abcorr
*@param   obsrvr               Name of the observing body. 
*@param   step                 Step size in seconds for finding FOV events. 
*@param   cnfine               SPICE window to which the search is restricted. 
*cfntime is a spiceDouble window. It must be initalized with wninsd()
*@param	  etstart optional... can be used in place of a spicewindow
*@param   etstop

*@returns A 2-demesional ruby array of windows
*@returns [[start,end],[start,end]...]
*@returns times are returned in et
*
*/
VALUE gftfov(int argc, VALUE *argv, VALUE self) {
	SpiceInt i, niv;
	double beg, end;
	SpiceCell *ptr;
	VALUE result;
	SPICEDOUBLE_CELL(intv, 5000);
	SPICEDOUBLE_CELL(win, 5000);

	if (argc < 8 || argc > 9){
		rb_raise(rb_eArgError, "need 8 parameters!");
		return Qnil;
	}

	Check_Type(argv[0], T_STRING);/*inst*/
	Check_Type(argv[1], T_STRING);/*target*/
	Check_Type(argv[2], T_STRING);/*tshape*/
	Check_Type(argv[3], T_STRING);/*tframe*/
	Check_Type(argv[4], T_STRING);/*abcorr*/
	Check_Type(argv[5], T_STRING);/*obsvr*/
	Check_Type(argv[6], T_FLOAT);/*step*/


	if (argc == 8){
	Check_Type(argv[7], T_DATA);/*cnfine*/
	Data_Get_Struct(argv[7], SpiceCell, ptr);

	gftfov_c (      StringValuePtr(argv[0]),
			StringValuePtr(argv[1]),
			StringValuePtr(argv[2]),
			StringValuePtr(argv[3]),
			StringValuePtr(argv[4]),
			StringValuePtr(argv[5]),
			NUM2DBL(argv[6]),
			ptr,
			&intv  );
	}
	else{
		Check_Type(argv[7], T_FLOAT);/*start et*/
		Check_Type(argv[8], T_FLOAT);/*stop et*/

		wninsd_c(NUM2DBL(argv[7]), NUM2DBL(argv[8]), &win);
		gftfov_c (      StringValuePtr(argv[0]),
				StringValuePtr(argv[1]),
				StringValuePtr(argv[2]),
				StringValuePtr(argv[3]),
				StringValuePtr(argv[4]),
				StringValuePtr(argv[5]),
				NUM2DBL(argv[6]),
				&win,
				&intv  );
	}

	check_spice_error();
	niv = wncard_c(&intv);
	if (niv == 0)
		return Qnil;

	result = rb_ary_new2(niv);

	for (i = 0; i < niv; i++){
		wnfetd_c(&intv, i, &beg, &end);
		rb_ary_push(result, rb_ary_new3(2, INT2FIX(beg), INT2FIX(end)));
	}

	return result;

}

/*!@fn VALUE gfrfov(int argc, VALUE *argv, VALUE self)
*@brief Determine time intervals when a specified ray intersects the 
*   space bounded by the field-of-view (FOV) of a specified 
*   instrument. 
*@brief http://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/gfrfov_c.html
*
*@param   inst                 Name of the instrument. 
*@param   raydir               Ray's direction vector. 
*@param   rframe               Reference frame of ray's direction vector. 
*@param   abcorr               Aberration correction flag. see @ref abcorr
*@param   obsrvr               Name of the observing body. 
*@param   step                 Step size in seconds for finding FOV events. 
*@param   cnfine               SPICE window to which the search is restricted. 
*cfntime is a spiceDouble window. It must be initalized with wninsd()
*
*@returns A 2-demesional ruby array of windows
*@returns [[start,end],[start,end]...]
*@returns times are returned in et
*
*/
VALUE gfrfov(int argc, VALUE *argv, VALUE self){ 
	SpiceInt i, niv;
	double beg, end, ray[3];
	VALUE result;
	SPICEDOUBLE_CELL(intv, 5000);
	SPICEDOUBLE_CELL(win, 5000);
	SpiceCell *ptr;

	if (argc < 7 || argc > 8){
		rb_raise(rb_eArgError, "need 7 parameters!");
		return Qnil;
	}

	Check_Type(argv[0], T_STRING);/*inst*/
	Check_Type(argv[1], T_ARRAY);/*raydir double[3]*/
	Check_Type(argv[2], T_STRING);/*rframe*/
	Check_Type(argv[3], T_STRING);/*abcorr*/
	Check_Type(argv[4], T_STRING);/*obsvr*/
	Check_Type(argv[5], T_FLOAT);/*step*/

	for (i=0; i < 3; i++) {
		Check_Type(RARRAY_PTR(argv[1])[i], T_FLOAT);
		ray[i] = NUM2DBL(RARRAY_PTR(argv[1])[i]);
	}
	
	if (argc == 7){
	Check_Type(argv[6], T_DATA);/*cnfine*/
	Data_Get_Struct(argv[6], SpiceCell, ptr);

	gfrfov_c(StringValuePtr(argv[0]),
			ray,
			StringValuePtr(argv[2]),
			StringValuePtr(argv[3]),
			StringValuePtr(argv[4]),
			NUM2DBL(argv[5]),
			ptr,
			&intv);
	}
	else{
		Check_Type(argv[6], T_FLOAT);/*start et*/
		Check_Type(argv[7], T_FLOAT);/*stop et*/

		wninsd_c(NUM2DBL(argv[6]), NUM2DBL(argv[7]), &win);

		gfrfov_c(StringValuePtr(argv[0]),
				ray,
				StringValuePtr(argv[2]),
				StringValuePtr(argv[3]),
				StringValuePtr(argv[4]),
				NUM2DBL(argv[5]),
				&win,
				&intv);


	}

	check_spice_error();
	niv = wncard_c(&intv);
	if (niv == 0)
		return Qnil;

	result = rb_ary_new2(niv);

	for (i = 0; i < niv; i++){
		wnfetd_c(&intv, i, &beg, &end);
		rb_ary_push(result, rb_ary_new3(2, INT2FIX(beg), INT2FIX(end)));
	}

	return result;

}
/*!@fn VALUE spkezp(int argc, VALUE *argv, VALUE self)
*@brief  Return the position of a target body relative to an observing 
*   body, optionally corrected for light time (planetary aberration) 
*   and stellar aberration.
*@brief http://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/spkezp_c.html
*
*@param   targ          Target body NAIF ID code. 
*@param   et            Observer epoch. 
*@param   ref           Reference frame of output position vector. 
*@param   abcorr        Aberration correction flag.  see @ref abcorr
*@param   obs           Observing body NAIF ID code. 
*
*@returns A ruby array containing the position vector, and light time between obs and targ
*@returns [[ptargX][ptargY][ptargZ], lt]
*/
VALUE spkezp(int argc, VALUE *argv, VALUE self) {
  double ptarg[3] = {0.0,0.0,0.0};
  double lt;
  VALUE rb_ptarg;
  SpiceInt i;

  if (argc != 5) {
    rb_raise(rb_eArgError, "need 5 parameters!");
    return Qnil;
  }

  Check_Type(argv[0], T_FIXNUM);
  Check_Type(argv[1], T_FLOAT);
  Check_Type(argv[2], T_STRING);
  Check_Type(argv[3], T_STRING);
  Check_Type(argv[4], T_FIXNUM);

  spkezp_c(FIX2INT(argv[0]), NUM2DBL(argv[1]), StringValuePtr(argv[2]), StringValuePtr(argv[3]), FIX2INT(argv[4]), ptarg, &lt);

  rb_ptarg = rb_ary_new();
  for (i=0; i < 3; i++)
    rb_ary_push(rb_ptarg, rb_float_new(ptarg[i]));

  check_spice_error();

  return rb_ary_new3(2, rb_ptarg, rb_float_new(lt));
}
/*!@fn VALUE spkezr(int argc, VALUE *argv, VALUE self)
*@brief  Return the state (position and velocity) of a target body 
*   relative to an observing body, optionally corrected for light 
*   time (planetary aberration) and stellar aberration.
*@brief http://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/spkezr_c.html
*
*@param   targ          Target body name. 
*@param   et            Observer epoch. 
*@param   ref           Reference frame of output state vector. 
*@param   abcorr        Aberration correction flag. See @ref abcorr 
*@param   obs           Observing body name. 
*
*@returns A ruby array, containing the State Vector, and light time between obs and targ
*@returns [[x,y,z,vx,vy,vz],lt]
*/
VALUE spkezr(int argc, VALUE *argv, VALUE self) {
  double starg[6] = {0.0,0.0,0.0,0.0,0.0,0.0};
  double lt;
  VALUE rb_starg;
  SpiceInt i;

  if (argc != 5) {
    rb_raise(rb_eArgError, "need 5 parameters!");
    return Qnil;
  }

  Check_Type(argv[0], T_STRING);
  Check_Type(argv[1], T_FLOAT);
  Check_Type(argv[2], T_STRING);
  Check_Type(argv[3], T_STRING);
  Check_Type(argv[4], T_STRING);

  spkezr_c(StringValuePtr(argv[0]), NUM2DBL(argv[1]), StringValuePtr(argv[2]), StringValuePtr(argv[3]), StringValuePtr(argv[4]), starg, &lt);

  rb_starg = rb_ary_new();
  for (i=0; i < 6; i++)
    rb_ary_push(rb_starg, rb_float_new(starg[i]));

  check_spice_error();

  return rb_ary_new3(2, rb_starg, rb_float_new(lt));
}
/*!@fn VALUE vdist(VALUE self, VALUE v1, VALUE v2)
*@brief   Return the distance between two three-dimensional vectors.
*@brief http://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/vdist_c.html
*
*@param v1 
*passed as a 3 item ruby array
*@param v2
*passed as a 3 item ruby array
*
*@returns A float representing the distance
*/
VALUE vdist(VALUE self, VALUE v1, VALUE v2) {
  double vec1[3], vec2[3], res;
  SpiceInt i;

  Check_Type(v1, T_ARRAY);
  Check_Type(v2, T_ARRAY);

  if (RARRAY_LEN(v1) != 3) {
    rb_raise(rb_eArgError, "The array should have 3 items in it");
    return Qnil;
  }

  if (RARRAY_LEN(v2) != 3) {
    rb_raise(rb_eArgError, "The array should have 3 items in it");
    return Qnil;
  }

  for (i=0; i < 3; i++) {
    Check_Type(RARRAY_PTR(v1)[i], T_FLOAT);
    vec1[i] = NUM2DBL(RARRAY_PTR(v1)[i]);
  }

  for (i=0; i < 3; i++) {
    Check_Type(RARRAY_PTR(v2)[i], T_FLOAT);
    vec2[i] = NUM2DBL(RARRAY_PTR(v2)[i]);
  }

  res = vdist_c(vec1, vec2);

  check_spice_error();

  return rb_float_new(res);
}
/*!@fn VALUE bodvcd(VALUE self, VALUE bodyid, VALUE item)
*@brief Fetch from the kernel pool the double precision values of an item
*   associated with a body, where the body is specified by an integer ID
*   code.
*@brief http://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/bodvcd_c.html
*
*@param   bodyid        Body ID code.
*@param   item          Item for which values are desired. ("RADII","NUT_PREC_ANGLES", etc. ) 
*
*@returns a ruby array containing the values
*
*/
VALUE bodvcd(VALUE self, VALUE bodyid, VALUE item) {
  double values[32];
  SpiceInt vret, i;
  VALUE rb_values;

  Check_Type(bodyid, T_FIXNUM);
  Check_Type(item, T_STRING);

  bodvcd_c(NUM2INT(bodyid), StringValuePtr(item), 32, &vret, values);

  rb_values = rb_ary_new2(vret);
  for (i=0; i<vret; i++)
    rb_ary_push(rb_values, rb_float_new(values[i]));

  check_spice_error();

  return rb_values;
}
/*!@fn VALUE gdpool(VALUE self, VALUE name, VALUE start)
*@brief  Return the d.p. value of a kernel variable from the kernel pool.
*@brief http://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/gdpool_c.html
*
*@param   name          Name of the variable whose value is to be returned.
*@param   start         Which component to start retrieving for name
*
*@returns nil if the varible was not found
*@returns A ruby array of values otherwise
*/
VALUE gdpool(VALUE self, VALUE name, VALUE start) {
  double values[32];
  SpiceInt vret, i;
  VALUE rb_values;
  SpiceBoolean found;

  Check_Type(name, T_STRING);
  Check_Type(start, T_FIXNUM);

  gdpool_c(StringValuePtr(name), NUM2INT(start), 32, &vret, values, &found);

  if (found == 0) {
    rb_raise(rb_eRuntimeError, "Variable not found in the pool");
    return Qnil;
  }

  rb_values = rb_ary_new2(vret);
  for (i=0; i<vret; i++)
    rb_ary_push(rb_values, rb_float_new(values[i]));

  check_spice_error();

  return rb_values;
}
/*!@fn VALUE lspcn(VALUE self, VALUE body, VALUE et, VALUE abcorr)
*@brief   Compute L_s, the planetocentric longitude of the sun, as seen 
*   from a specified body.
*@brief http://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/lspcn_c.html
* 
*@param   body          Name of central body. 
*@param   et            Epoch in seconds past J2000 TDB. 
*@param   abcorr        Aberration correction. See @ref abcorr
*
*/
VALUE lspcn(VALUE self, VALUE body, VALUE et, VALUE abcorr) {
  SpiceDouble res;

  Check_Type(body, T_STRING);
  Check_Type(et, T_FLOAT);
  Check_Type(abcorr, T_STRING);

  res = lspcn_c(StringValuePtr(body), NUM2DBL(et), StringValuePtr(abcorr));

  check_spice_error();

  return rb_float_new(res);
}
/*!@fn VALUE rpd(VALUE self)
*@brief    Return the number of radians per degree.
*@brief http://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/rpd_c.html
*@param NONE
*
*@returns acos(-1.0)/180
*/
VALUE rpd(VALUE self) {
  return rb_float_new(rpd_c());
}
/*!@fn VALUE latrec(VALUE self, VALUE radius, VALUE lon, VALUE lat)
*@brief  Convert from latitudinal coordinates to rectangular coordinates.
*@brief http://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/latrec_c.html
*
*@param   radius        Distance of a point from the origin.
*@param   longitude     Longitude of point in radians.
*@param   latitude      Latitude of point in radians.
*
*@returns A ruby array containing the rectangular vector. EX: [X,Y,Z]
*
*/
VALUE latrec(VALUE self, VALUE radius, VALUE lon, VALUE lat) {
  double spoint[3];

  Check_Type(radius, T_FLOAT);
  Check_Type(lon, T_FLOAT);
  Check_Type(lat, T_FLOAT);

  latrec_c(NUM2DBL(radius), NUM2DBL(lon), NUM2DBL(lat), spoint);

  check_spice_error();

  return rb_ary_new3(3, rb_float_new(spoint[0]), rb_float_new(spoint[1]), rb_float_new(spoint[2]));
}
/*!@fn VALUE vsub(VALUE self, VALUE rb_v1, VALUE rb_v2)
*@brief Compute the difference between two 3-dimensional, double 
*   precision vectors. 
*@brief http://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/vsub_c.html
*
*@paran v1 First vector (minuend).\n
*passed in as a 3 item ruby array
*@param v2 Second vector (subtrahend).\n
*passed in as a 3 item ruby array
*
*@returns A 3 item ruby array representing the difference vector
*
*/
VALUE vsub(VALUE self, VALUE rb_v1, VALUE rb_v2) {
  double v1[3], v2[3], vd[3];
  SpiceInt i;

  Check_Type(rb_v1, T_ARRAY);
  Check_Type(rb_v2, T_ARRAY);
  
  for(i = 0; i < 3; i++) {
    Check_Type(RARRAY_PTR(rb_v1)[i], T_FLOAT);
    Check_Type(RARRAY_PTR(rb_v2)[i], T_FLOAT);
    v1[i] = NUM2DBL(RARRAY_PTR(rb_v1)[i]);
    v2[i] = NUM2DBL(RARRAY_PTR(rb_v2)[i]);
  }

  vsub_c(v1, v2, vd);

  check_spice_error();

  return rb_ary_new3(3, rb_float_new(vd[0]), rb_float_new(vd[1]), rb_float_new(vd[2]));
}
/*!@fn VALUE ucrss(VALUE self, VALUE rb_v1, VALUE rb_v2)
* @brief   Compute the normalized cross product of two 3-vectors. 
* @brief http://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/ucrss_c.html
*
*@param v1              Left vector for cross product. 
*@param v2              Right vector for cross product. 
*
*@returns a 3 item ruby array representing the Normalized cross product (v1xv2) / |v1xv2|. 
*
*/
VALUE ucrss(VALUE self, VALUE rb_v1, VALUE rb_v2) {
  double v1[3], v2[3], vd[3];
  SpiceInt i;

  Check_Type(rb_v1, T_ARRAY);
  Check_Type(rb_v2, T_ARRAY);
  
  for(i = 0; i < 3; i++) {
    Check_Type(RARRAY_PTR(rb_v1)[i], T_FLOAT);
    Check_Type(RARRAY_PTR(rb_v2)[i], T_FLOAT);
    v1[i] = NUM2DBL(RARRAY_PTR(rb_v1)[i]);
    v2[i] = NUM2DBL(RARRAY_PTR(rb_v2)[i]);
  }

  ucrss_c(v1, v2, vd);

  check_spice_error();

  return rb_ary_new3(3, rb_float_new(vd[0]), rb_float_new(vd[1]), rb_float_new(vd[2]));
}
/*!@fn VALUE vcrss(VALUE self, VALUE rb_v1, VALUE rb_v2)
*@brief Compute the cross product of two 3-dimensional vectors.
*@brief http://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/vcrss_c.html
*
*@param v1  Left hand vector for cross product.
*@param v2  Right hand vector for cross product.
*
*@returns a 3 item ruby array representing the Cross product v1xv2.
*
*/
VALUE vcrss(VALUE self, VALUE rb_v1, VALUE rb_v2) {
  double v1[3], v2[3], vd[3];
  SpiceInt i;

  Check_Type(rb_v1, T_ARRAY);
  Check_Type(rb_v2, T_ARRAY);
  
  for(i = 0; i < 3; i++) {
    Check_Type(RARRAY_PTR(rb_v1)[i], T_FLOAT);
    Check_Type(RARRAY_PTR(rb_v2)[i], T_FLOAT);
    v1[i] = NUM2DBL(RARRAY_PTR(rb_v1)[i]);
    v2[i] = NUM2DBL(RARRAY_PTR(rb_v2)[i]);
  }

  vcrss_c(v1, v2, vd);

  check_spice_error();

  return rb_ary_new3(3, rb_float_new(vd[0]), rb_float_new(vd[1]), rb_float_new(vd[2]));
}
/*!@fn VALUE vdot(VALUE self, VALUE rb_v1, VALUE rb_v2)
*@brief Compute the dot product of two double precision, 3-dimensional vectors.
*@brief http://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/vdot_c.html
*
*@param v1 First vector in the dot product.
*@param v2 Second vector in the dot product.
*
*@returns A float representing the value of the dot product of v1 and v2.
*
*/
VALUE vdot(VALUE self, VALUE rb_v1, VALUE rb_v2) {
  double v1[3], v2[3], res;
  SpiceInt i;

  Check_Type(rb_v1, T_ARRAY);
  Check_Type(rb_v2, T_ARRAY);
  
  for(i = 0; i < 3; i++) {
    Check_Type(RARRAY_PTR(rb_v1)[i], T_FLOAT);
    Check_Type(RARRAY_PTR(rb_v2)[i], T_FLOAT);
    v1[i] = NUM2DBL(RARRAY_PTR(rb_v1)[i]);
    v2[i] = NUM2DBL(RARRAY_PTR(rb_v2)[i]);
  }

  res = vdot_c(v1, v2);

  check_spice_error();

  return rb_float_new(res);
}
/*!@fn VALUE sctiks(VALUE self, VALUE sc, VALUE slkstr)
*@brief  Convert a spacecraft clock format string to number of "ticks". 
*@http://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/sctiks_c.html
*
*@param sc            NAIF spacecraft identification code. 
*@param clkstr        Character representation of a spacecraft clock. 
*
*@returns A float representing the numbe of ticks
*
*/
VALUE sctiks(VALUE self, VALUE sc, VALUE slkstr) {
  double ticks;

  Check_Type(sc, T_FIXNUM);
  Check_Type(slkstr, T_STRING);

  sctiks_c(NUM2INT(sc), StringValuePtr(slkstr), &ticks);

  check_spice_error();

  return rb_float_new(ticks);
}
/*!@fn VALUE vsep(VALUE self, VALUE rb_v1, VALUE rb_v2)
*@brief Find the separation angle in radians between two double 
*   precision, 3-dimensional vectors.  This angle is defined as zero 
*   if either vector is zero. 
*@brief http://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/vsep_c.html
*
*@param v1    first vector
*passed in as a 3 item ruby array
*@param v2    second vector
*passed in as a 3 item ruby array
*
*@returns A float representing the seperation angle in radians
*/
VALUE vsep(VALUE self, VALUE rb_v1, VALUE rb_v2) {
  double v1[3], v2[3], res;
  SpiceInt i;

  Check_Type(rb_v1, T_ARRAY);
  Check_Type(rb_v2, T_ARRAY);
  
  for(i = 0; i < 3; i++) {
    Check_Type(RARRAY_PTR(rb_v1)[i], T_FLOAT);
    Check_Type(RARRAY_PTR(rb_v2)[i], T_FLOAT);
    v1[i] = NUM2DBL(RARRAY_PTR(rb_v1)[i]);
    v2[i] = NUM2DBL(RARRAY_PTR(rb_v2)[i]);
  }

  res = vsep_c(v1, v2);

  check_spice_error();

  return rb_float_new(res);
}
/*!@fn VALUE ckgp(VALUE self, VALUE inst, VALUE sclkdp, VALUE tol, VALUE ref)
*@brief Get pointing (attitude) for a specified spacecraft clock time.
*@brief http://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/ckgp_c.html
*
*@param   inst          NAIF ID of instrument, spacecraft, or structure.
*@param   sclkdp        Encoded spacecraft clock time. 
*@param   tol           Time tolerance. 
*@param   ref           Reference frame. 
*
*@returns 'nil' if no pointing information was found, otherwise:
*
*@returns A ruby array containing the c-matrix as it's first element and the sclk time as the second
*
*For Reference:
*@code return rb_ary_new3(2, rb_ary_new3(3, rb_ary_new3(3, rb_float_new(cmat[0][0]), rb_float_new(cmat[0][1]), rb_float_new(cmat[0][2])), rb_ary_new3(3, rb_float_new(cmat[1][0]), rb_float_new(cmat[1][1]), rb_float_new(cmat[1][2])), rb_ary_new3(3, rb_float_new(cmat[2][0]), rb_float_new(cmat[2][1]), rb_float_new(cmat[2][2]))), rb_float_new(slkout));@endcode
*
*/
VALUE ckgp(VALUE self, VALUE inst, VALUE sclkdp, VALUE tol, VALUE ref) {
  double cmat[3][3];
  double slkout;
  SpiceBoolean found;

  Check_Type(inst, T_FIXNUM);
  Check_Type(sclkdp, T_FLOAT);
  Check_Type(tol, T_FLOAT);
  Check_Type(ref, T_STRING);

  ckgp_c(NUM2INT(inst), NUM2DBL(sclkdp), NUM2DBL(tol), StringValuePtr(ref), cmat, &slkout, &found);

  if (found == SPICEFALSE) {
    rb_raise(rb_eRuntimeError, "Pointing info not found for this sc time");
    return Qnil;
  }

  check_spice_error();

  return rb_ary_new3(2, rb_ary_new3(3, rb_ary_new3(3, rb_float_new(cmat[0][0]), rb_float_new(cmat[0][1]), rb_float_new(cmat[0][2])), rb_ary_new3(3, rb_float_new(cmat[1][0]), rb_float_new(cmat[1][1]), rb_float_new(cmat[1][2])), rb_ary_new3(3, rb_float_new(cmat[2][0]), rb_float_new(cmat[2][1]), rb_float_new(cmat[2][2]))), rb_float_new(slkout));
}
/*!@fn VALUE mxv(VALUE self, VALUE rb_m1, VALUE rb_v1)
*@brief Multiply a 3x3 double precision matrix with a 3-dimensional
*   double precision vector.
*@brief http://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/mxv_c.html
*
*@param m1 3x3 double precision matrix.
* passed in as ruby array [3][3]
*@param v1 3-dimensional double precision vector.
* passed in as a 3 item ruby array
*
*@returns A 3 item ruby array representing the product m1*vin.
*
*/
VALUE mxv(VALUE self, VALUE rb_m1, VALUE rb_v1) {
  double v2[3];
  double v1[3];
  double m1[3][3];
  SpiceInt i, j;

  if (RARRAY_LEN(rb_v1) != 3) {
    rb_raise(rb_eArgError, "V1 should have 3 items in it");
    return Qnil;
  }

  if (RARRAY_LEN(rb_m1) != 3) {
    rb_raise(rb_eArgError, "M1 should be 3x3");
    return Qnil;
  }

  for(i = 0; i < 3; i++) {
    Check_Type(RARRAY_PTR(rb_v1)[i], T_FLOAT);
    v1[i] = NUM2DBL(RARRAY_PTR(rb_v1)[i]);

    if (RARRAY_LEN(RARRAY_PTR(rb_m1)[i]) != 3) {
      rb_raise(rb_eArgError, "M1 should be 3x3.");
      return Qnil;
    }

    for(j = 0; j < 3; j++) {
      Check_Type(RARRAY_PTR(RARRAY_PTR(rb_m1)[i])[j], T_FLOAT);
      m1[i][j] = NUM2DBL(RARRAY_PTR(RARRAY_PTR(rb_m1)[i])[j]);
    }
  }

  mxv_c(m1, v1, v2);

  check_spice_error();

  return rb_ary_new3(3, rb_float_new(v2[0]), rb_float_new(v2[1]), rb_float_new(v2[2]));
}

/*!@fn VALUE gfoclt(int argc, VALUE* argv, VALUE self)
*@brief Determine time intervals when an observer sees one target 
*   occulted by, or in transit across, another.
*@brief http://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/gfoclt_c.html
*
*@param   occtyp              Type of occultation.\n
*Supported values:\n
*"FULL"\n 
*"ANNULAR"\n   
*"PARTIAL"\n      
*"ANY"\n\n
*@param   front               Name of body occulting the other. 
*@param   fshape              Type of shape model used for front body. \n
*Supported values:\n
*"ELLIPSOID"\n
* "POINT"\n\n
*@param   fframe              Body-fixed, body-centered frame for front body. 
*@param   back                Name of body occulted by the other. 
*@param   bshape              Type of shape model used for back body.\n 
*Supported values:\n
*"ELLIPSOID"\n
* "POINT"\n\n
*@param   bframe              Body-fixed, body-centered frame for back body. 
*@param   abcorr              Aberration correction flag. see @ref abcorr 
*@param   obsrvr              Name of the observing body. 
*@param   step                Step size in seconds for finding occultation events. 
*@param   cnfine              SPICE window to which the search is restricted. 
*cfntime is a spiceDouble window. It must be initalized with wninsd()
*
*@returns A 2-demesional ruby array of windows
*@returns [[start,end],[start,end]...]
*@returns times are returned in et
*/
VALUE gfoclt(int argc, VALUE* argv, VALUE self) {					
	SPICEDOUBLE_CELL (result, 5000);
	SpiceInt i, niv;
	double beg, end;
	VALUE rb_result;
	SpiceCell *cnfine;
	SPICEDOUBLE_CELL(win, 5000);

	if (argc < 11 || argc > 12){
		rb_raise(rb_eArgError, "need 11 parameters!");
		return Qnil;
	}
	Check_Type(argv[0], T_STRING); /*occtyp*/
	Check_Type(argv[1], T_STRING); /*front*/
	Check_Type(argv[2], T_STRING); /*fshape*/
	Check_Type(argv[3], T_STRING); /*fframe*/
	Check_Type(argv[4], T_STRING); /*back*/
	Check_Type(argv[5], T_STRING); /*bshape*/
	Check_Type(argv[6], T_STRING); /*bframe*/
	Check_Type(argv[7], T_STRING); /*abcorr*/
	Check_Type(argv[8], T_STRING); /*obsrvr*/
	Check_Type(argv[9], T_FLOAT); /*step*/


	if (argc == 11){
		Check_Type(argv[10],T_DATA); /*cnfine*/
		Data_Get_Struct(argv[10], SpiceCell, cnfine);

		gfoclt_c( StringValuePtr(argv[0]), 
				StringValuePtr(argv[1]),  
				StringValuePtr(argv[2]),
				StringValuePtr(argv[3]), 
				StringValuePtr(argv[4]), 
				StringValuePtr(argv[5]), 
				StringValuePtr(argv[6]), 
				StringValuePtr(argv[7]),  
				StringValuePtr(argv[8]),  
				NUM2DBL(argv[9]), 
				cnfine,
				&result);
	}
	else{
		Check_Type(argv[10], T_FLOAT);/*start et*/
		Check_Type(argv[11], T_FLOAT);/*stop et*/

		wninsd_c(NUM2DBL(argv[10]), NUM2DBL(argv[11]), &win);
		
		gfoclt_c( StringValuePtr(argv[0]), 
				StringValuePtr(argv[1]),  
				StringValuePtr(argv[2]),
				StringValuePtr(argv[3]), 
				StringValuePtr(argv[4]), 
				StringValuePtr(argv[5]), 
				StringValuePtr(argv[6]), 
				StringValuePtr(argv[7]),  
				StringValuePtr(argv[8]),  
				NUM2DBL(argv[9]), 
				&win,
				&result);
	}

	check_spice_error();

	niv = wncard_c(&result);
	if (niv == 0)
		return Qnil;

	rb_result = rb_ary_new2(niv);

	for (i = 0; i < niv; i++) {
		wnfetd_c(&result, i, &beg, &end);
		rb_ary_push(rb_result, rb_ary_new3(2, INT2FIX(beg), INT2FIX(end)));
	}

	return rb_result;
}
/*!@fn VALUE xf2eul(VALUE self, VALUE rb_cmat, VALUE rb_axisa, VALUE rb_axisb, VALUE rb_axisc)
*@brief Convert a state transformation matrix to Euler angles and their 
   derivatives with respect to a specified set of axes. 
*@brief http://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/xf2eul_c.html 
*
*@param   xform         A state transformation matrix. 
*passed in by a ruby array [6][6]
*@param   axisa         Axis A of the Euler angle factorization. 
*@param   axisb         Axis B of the Euler angle factorization. 
*@param   axisc         Axis C of the Euler angle factorization. 
*
*@returns A ruby array contatining (An array of Euler angles and their derivatives) and wether or not it is a unique representation(true|false) 
*@returns [[6],unique]
*@returns For Reference:
*@code   return rb_ary_new3(2, rb_ary_new3(6, rb_float_new(eulang[0]), rb_float_new(eulang[1]), rb_float_new(eulang[2]), rb_float_new(eulang[3]), rb_float_new(eulang[4]), rb_float_new(eulang[5])), rb_unique);@endcode
*/
VALUE xf2eul(VALUE self, VALUE rb_cmat, VALUE rb_axisa, VALUE rb_axisb, VALUE rb_axisc) {
  double cmat[6][6];
  double eulang[6];
  SpiceBoolean unique;
  SpiceInt i, j;
  VALUE rb_unique;

  if (RARRAY_LEN(rb_cmat) != 6) {
    rb_raise(rb_eArgError, "cmat should be 6x6");
    return Qnil;
  }

  for(i = 0; i < 6; i++) {
    if (RARRAY_LEN(RARRAY_PTR(rb_cmat)[i]) != 6) {
      rb_raise(rb_eArgError, "cmat should be 6x6.");
      return Qnil;
    }

    for(j = 0; j < 6; j++) {
      Check_Type(RARRAY_PTR(RARRAY_PTR(rb_cmat)[i])[j], T_FLOAT);
      cmat[i][j] = NUM2DBL(RARRAY_PTR(RARRAY_PTR(rb_cmat)[i])[j]);
    }
  }

  xf2eul_c(cmat, NUM2INT(rb_axisa), NUM2INT(rb_axisb), NUM2INT(rb_axisc), eulang, &unique);

  check_spice_error();

  if (unique == SPICETRUE)
    rb_unique = Qtrue;
  else
    rb_unique = Qfalse;

  return rb_ary_new3(2, rb_ary_new3(6, rb_float_new(eulang[0]), rb_float_new(eulang[1]), rb_float_new(eulang[2]), rb_float_new(eulang[3]), rb_float_new(eulang[4]), rb_float_new(eulang[5])), rb_unique);
}
/*!@fn VALUE m2eul(VALUE self, VALUE rb_r, VALUE rb_axis3, VALUE rb_axis2, VALUE rb_axis1)
*@brief Factor a rotation matrix as a product of three rotations about
*   specified coordinate axes.
*@brief http://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/m2eul_c.html
*
*@param r A rotation matrix to be factored.
*passed in as a ruby array [3][3]
*@param axis1 Numbers of third, second, and first rotation axes.
*@param axis2
*@param axis3
*
*@returns A ruby array representing the 3 Euler angles.
*@returns [angle3, angle2, angle1]
*
*/
VALUE m2eul(VALUE self, VALUE rb_r, VALUE rb_axis3, VALUE rb_axis2, VALUE rb_axis1) {
  double r[3][3];
  SpiceDouble angle3, angle2, angle1;
  SpiceInt i, j;

  if (RARRAY_LEN(rb_r) != 3) {
    rb_raise(rb_eArgError, "R should be 3x3");
    return Qnil;
  }

  for(i=0; i<3; i++) {
    if (RARRAY_LEN(RARRAY_PTR(rb_r)[i]) !=3) {
      rb_raise(rb_eArgError, "R should be 3x3");
      return Qnil;
    }

    for (j=0; j<3; j++) {
      Check_Type(RARRAY_PTR(RARRAY_PTR(rb_r)[i])[j], T_FLOAT);
      r[i][j] = NUM2DBL(RARRAY_PTR(RARRAY_PTR(rb_r)[i])[j]);
    }
  }

  m2eul_c(r, NUM2INT(rb_axis3), NUM2INT(rb_axis2), NUM2INT(rb_axis1), &angle3, &angle2, &angle1);

  check_spice_error();

  return rb_ary_new3(3, rb_float_new(angle3), rb_float_new(angle2), rb_float_new(angle1));
}

/*!@fn VALUE m2q(VALUE self, VALUE from, VALUE to, VALUE et)
*@brief    Find the rotation matrix corresponding to a specified unit quaternion. 
*@brief http://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/m2q_c.html
*
*@param   r          A rotation matrix (In ruby [[3],[3],[3]])
*
*@returns a ruby array[4] containing A unit quaternion representing `r'.
*
*/
VALUE m2q(int argc, VALUE *argv, VALUE self){
  double r[3][3];
  double q[4];
  int i, j;

  if (argc != 1) {
    rb_raise(rb_eArgError, "Need 1 parameter!");
    return Qnil;
  }
  Check_Type(argv[0], T_ARRAY);

  if (RARRAY_LEN(argv[0]) != 3) {
    rb_raise(rb_eArgError, "r should be 3x3");
    return Qnil;
  }

  for(i=0; i<3; i++) {
    if (RARRAY_LEN(RARRAY_PTR(argv[0])[i]) !=3) {
      rb_raise(rb_eArgError, "r should be 3x3");
      return Qnil;
    }
    for (j=0; j<3; j++) {
      Check_Type(RARRAY_PTR(RARRAY_PTR(argv[0])[i])[j], T_FLOAT);
      r[i][j] = NUM2DBL(RARRAY_PTR(RARRAY_PTR(argv[0])[i])[j]);
    }
  }

  m2q_c(r, q);

  check_spice_error();

  return rb_ary_new3(4, rb_float_new(q[0]),  rb_float_new(q[1]), rb_float_new(q[2]), rb_float_new(q[3]));
}

/*!@fn VALUE sce2c(VALUE self, VALUE sc, VALUE et)
*@brief Convert ephemeris seconds past J2000 (ET) to continuous encoded  
*   spacecraft clock (`ticks').  Non-integral tick values may be returned. 
*@brief http://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/sce2c_c.html
*
*@param  sc         NAIF spacecraft ID code.
*@param  et         Ephemeris time, seconds past J2000.
*
*@returns A float representing: SCLK, encoded as ticks since spacecraft clock start.
*/
VALUE sce2c(VALUE self, VALUE sc, VALUE et) {
  double sctk;

  Check_Type(sc, T_FIXNUM);
  Check_Type(et, T_FLOAT);

  sce2c_c(NUM2INT(sc), NUM2DBL(et), &sctk);

  check_spice_error();

  return rb_float_new(sctk);
}
/*!@fn VALUE scs2e(VALUE self, VALUE sc, VALUE sclk)
*@brief Convert a spacecraft clock string to ephemeris seconds past J2000 (ET). 
*@brief http://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/scs2e_c.html
*
*@param sc           NAIF integer code for a spacecraft. 
*@param sclkch       An SCLK string.
*
*@returns A float representing; et-Ephemeris time, seconds past J2000. 
*
*/
VALUE scs2e(VALUE self, VALUE sc, VALUE sclk) {
  double et;

  Check_Type(sc, T_FIXNUM);
  Check_Type(sclk, T_STRING);

  scs2e_c(NUM2INT(sc), StringValuePtr(sclk), &et);

  check_spice_error();

  return rb_float_new(et);
}
/*!@fn VALUE sxform(VALUE self, VALUE from, VALUE to, VALUE et)
*@brief Return the state transformation matrix from one frame to 
*   another at a specified epoch. 
*@brief http://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/sxform_c.html
*
*@param   from          Name of the frame to transform from.
*@param   to            Name of the frame to transform to.
*@param   et            Epoch of the state transformation matrix.
*
*@returns A ruby array[6][6] representing the state transformation matrix.
*/
VALUE sxform(VALUE self, VALUE from, VALUE to, VALUE et) {
  double cmat[6][6];

  Check_Type(from, T_STRING);
  Check_Type(to, T_STRING);
  Check_Type(et, T_FLOAT);

  sxform_c(StringValuePtr(from), StringValuePtr(to), NUM2DBL(et), cmat);

  check_spice_error();

  return rb_ary_new3(6, rb_ary_new3(6, rb_float_new(cmat[0][0]), rb_float_new(cmat[0][1]), rb_float_new(cmat[0][2]), rb_float_new(cmat[0][3]), rb_float_new(cmat[0][4]), rb_float_new(cmat[0][5])),
			rb_ary_new3(6, rb_float_new(cmat[1][0]), rb_float_new(cmat[1][1]), rb_float_new(cmat[1][2]), rb_float_new(cmat[1][3]), rb_float_new(cmat[1][4]), rb_float_new(cmat[1][5])),
			rb_ary_new3(6, rb_float_new(cmat[2][0]), rb_float_new(cmat[2][1]), rb_float_new(cmat[2][2]), rb_float_new(cmat[2][3]), rb_float_new(cmat[2][4]), rb_float_new(cmat[2][5])),
			rb_ary_new3(6, rb_float_new(cmat[3][0]), rb_float_new(cmat[3][1]), rb_float_new(cmat[3][2]), rb_float_new(cmat[3][3]), rb_float_new(cmat[3][4]), rb_float_new(cmat[3][5])),
			rb_ary_new3(6, rb_float_new(cmat[4][0]), rb_float_new(cmat[4][1]), rb_float_new(cmat[4][2]), rb_float_new(cmat[4][3]), rb_float_new(cmat[4][4]), rb_float_new(cmat[4][5])),
			rb_ary_new3(6, rb_float_new(cmat[5][0]), rb_float_new(cmat[5][1]), rb_float_new(cmat[5][2]), rb_float_new(cmat[5][3]), rb_float_new(cmat[5][4]), rb_float_new(cmat[5][5])));
}
/*!@fn VALUE pxform(VALUE self, VALUE from, VALUE to, VALUE et)
*@brief   Return the matrix that transforms position vectors from one 
*   specified frame to another at a specified epoch.
*@brief http://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/pxform_c.html
*
*@param   from          Name of the frame to transform from.
*@param   to            Name of the frame to transform to.
*@param   et            Epoch of the rotation matrix.
*
*@returns a ruby array[3][3] containing the rotation matrix
*
*/
VALUE pxform(VALUE self, VALUE from, VALUE to, VALUE et) {
  double cmat[3][3];

  Check_Type(from, T_STRING);
  Check_Type(to, T_STRING);
  Check_Type(et, T_FLOAT);

  pxform_c(StringValuePtr(from), StringValuePtr(to), NUM2DBL(et), cmat);

  check_spice_error();

  return rb_ary_new3(3, rb_ary_new3(3, rb_float_new(cmat[0][0]), rb_float_new(cmat[0][1]), rb_float_new(cmat[0][2])), rb_ary_new3(3, rb_float_new(cmat[1][0]), rb_float_new(cmat[1][1]), rb_float_new(cmat[1][2])), rb_ary_new3(3, rb_float_new(cmat[2][0]), rb_float_new(cmat[2][1]), rb_float_new(cmat[2][2])));
}

/*!@fn VALUE q2m(VALUE self, VALUE from, VALUE to, VALUE et)
*@brief    Find the rotation matrix corresponding to a specified unit quaternion. 
*@brief http://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/q2m_c.html
*
*@param   q          A unit quaternion. (in ruby a 4 element array)
*
*@returns a ruby array[3][3] containing the rotation matrix corresponding to 'q'
*
*/
VALUE q2m(int argc, VALUE *argv, VALUE self){
  double q[4];
  double r[3][3];
  int i;

  if (argc != 1) {
    rb_raise(rb_eArgError, "Need 1 parameter!");
    return Qnil;
  }
  Check_Type(argv[0], T_ARRAY);

  for (i=0; i < 4; i++){
    Check_Type(RARRAY_PTR(argv[0])[i], T_FLOAT);
    q[i] = NUM2DBL(RARRAY_PTR(argv[0])[i]);
  }

  q2m_c(q, r);

  check_spice_error();


  return rb_ary_new3(3, rb_ary_new3(3, rb_float_new(r[0][0]),  rb_float_new(r[0][1]), rb_float_new(r[0][2])), rb_ary_new3(3, rb_float_new(r[1][0]),  rb_float_new(r[1][1]), rb_float_new(r[1][2])), rb_ary_new3(3, rb_float_new(r[2][0]),  rb_float_new(r[2][1]), rb_float_new(r[2][2])));
}

/*!@fn VALUE recrad(int argc, VALUE* argv, VALUE self)
*@brief Convert rectangular coordinates to range, right ascension, and declination.
*@brief http://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/recrad_c.html
*
*@param rectan 	Rectangular coordinates of a point.
*passed in as a 3 item ruby array
*
*@returns A 3 item ruby array, representing: [range, ra. dec]
* @returns range         Distance of the point from the origin. 
* @returns ra            Right ascension in radians. 
* @return  dec           Declination in radians. 
*
*/
VALUE recrad(int argc, VALUE* argv, VALUE self) {
  double spoint[3];
  SpiceInt i;
  double range, ra, dec;

  if (argc != 1) {
    rb_raise(rb_eArgError, "Need just 1 parameter");
    return Qnil;
  }

  Check_Type(argv[0], T_ARRAY);

  if (RARRAY_LEN(argv[0]) != 3) {
    rb_raise(rb_eArgError, "The array should have 3 items in it");
    return Qnil;
  }

  for (i=0; i < 3; i++) {
    Check_Type(RARRAY_PTR(argv[0])[i], T_FLOAT);
    spoint[i] = NUM2DBL(RARRAY_PTR(argv[0])[i]);
  }

  recrad_c(spoint, &range, &ra, &dec);

  check_spice_error();

  return rb_ary_new3(3, rb_float_new(range), rb_float_new(ra), rb_float_new(dec));
}
/*!@fn VALUE spkpos(VALUE self, VALUE rb_targ, VALUE rb_et, VALUE rb_ref, VALUE rb_abcorr, VALUE rb_obs)
*@brief    Return the position of a target body relative to an observing 
*   body, optionally corrected for light time (planetary aberration) 
*   and stellar aberration.
*@brief http://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/spkpos_c.html
*
*@param   targ          Target body name. 
*@param   et            Observer epoch. 
*@param   ref           Reference frame of output position vector. 
*@param   abcorr        Aberration correction flag. see @ref abcorr 
*@param   obs           Observing body name.
*
*@returns A ruby array containing the target position vector and One way light time between observer and target.
*@returns [[X,Y,Z],lt]
*/
VALUE spkpos(VALUE self, VALUE rb_targ, VALUE rb_et, VALUE rb_ref, VALUE rb_abcorr, VALUE rb_obs) {
  SpiceDouble ptarg[3], lt;

  Check_Type(rb_targ, T_STRING);
  Check_Type(rb_et, T_FLOAT);
  Check_Type(rb_ref, T_STRING);
  Check_Type(rb_abcorr, T_STRING);
  Check_Type(rb_obs, T_STRING);

  spkpos_c(StringValuePtr(rb_targ), NUM2DBL(rb_et), StringValuePtr(rb_ref), StringValuePtr(rb_abcorr), StringValuePtr(rb_obs), ptarg, &lt);

  check_spice_error();

  return rb_ary_new3(2, rb_ary_new3(3, rb_float_new(ptarg[0]), rb_float_new(ptarg[1]), rb_float_new(ptarg[2])), rb_float_new(lt));
}

/*!@fn VALUE wninsd (int argc, VALUE *argv, VALUE self)
*@brief Insert an interval into a double precision window. 
*@brief http://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/wninsd_c.html
*
*@param   left       Left endpoint of new interval (et). 
*@param   right      Right endpoint of new interval (et). 
*
*@returns A pointer to the SpiceCell containing the window
*@returns This can then be passed into functions which require a spiceWindow
*/
#ifndef __APPLE__
VALUE wninsd (int argc, VALUE *argv, VALUE self) {
	VALUE rb_window;
	VALUE class = 0;
	SpiceCell *ptr = 0;
	SPICEDOUBLE_CELL (window, 5000);
	if (argc < 2 || argc > 3){
		rb_raise(rb_eArgError, "need 2 parameters!");
		return Qnil;
	}


	/*        !!!!!! Not ready to lose this bit incase we decide to rework 
	 *        !!!!!! Window handling at a later date to work on osX
	 *        !!!!!! This sets up a usable spiceDouble Cell without using
	 *        !!!!!! the NAIF marcos, which do some invisible magic
	 *        !!!!!! alicht@ser.asu.edu
	 *	//The extra 6's in here are defined in Spice_cel.h they 
	 *	//they are used to hold "spice cell control parameters
	 *	//In the case of DP cells, the only one set is [4], it's 
	 *	//set to the size of the cell minus the ctrl size (ie 5006 - 6 = 5000)
	 *        SpiceCell *window = malloc(sizeof(SpiceCell));
	 *        SpiceDouble windowArr[5006];
	 *        memset(windowArr, 0, 5006);
	 *
	 *        windowArr[4] = 5000;
	 *        window->dtype = SPICE_DP;
	 *        window->length = 0;
	 *        window->size = 5000;
	 *        window->card = SPICETRUE;
	 *        window->isSet = SPICEFALSE;
	 *        window->adjust = SPICEFALSE;
	 *        window->base = (void *)windowArr;
	 *        window->data = (void *)(&((SpiceDouble *) window->base)[6]);
	 */

	Check_Type(argv[0], T_FLOAT);
	Check_Type(argv[1], T_FLOAT);
	if ( argc == 3) {
          Check_Type(argv[2], T_DATA);
          Data_Get_Struct(argv[2], SpiceCell, ptr);
	  wninsd_c( NUM2DBL(argv[0]), NUM2DBL(argv[1]), ptr);
	}

	  wninsd_c( NUM2DBL(argv[0]), NUM2DBL(argv[1]), &window );

	check_spice_error();
        /*
	 *A spice cell is just a struct full of primitive types wrapped 
	 *around an array of the spice_type, an array of doubles in this 
	 *case
         */
	if (argc == 2)
	  rb_window = Data_Wrap_Struct(class, 0, 0, &window);
	else
	  rb_window = Data_Wrap_Struct(class, 0, 0, ptr);

	return rb_window;
}
#endif /*__APPLE__*/

VALUE wn2rb (VALUE self, VALUE window) {
	int niv, i;
	VALUE result;
	double beg, end;
	SpiceCell *ptr;

	Check_Type(window, T_DATA);
	Data_Get_Struct(window, SpiceCell, ptr);


	niv = wncard_c(ptr);
	if (niv == 0)
		return Qnil;

	result = rb_ary_new2(niv);

	for (i = 0; i < niv; i++){
		wnfetd_c(ptr, i, &beg, &end);
		rb_ary_push(result, rb_ary_new3(2, INT2FIX(beg), INT2FIX(end)));
	}

	return result;
}

VALUE mSpice;

void Init_spice(void) {
  mSpice = rb_define_module("Spice");
  
  rb_require("./spice_utils.rb"); 
  
  rb_define_module_function(mSpice, "bodn2c", bodn2c, -1);
  rb_define_module_function(mSpice, "bodvcd", bodvcd, 2);
  rb_define_module_function(mSpice, "ckgp", ckgp, 4);
  rb_define_module_function(mSpice, "dpr", dpr, -1);
  rb_define_module_function(mSpice, "eul2m", eul2m, -1);
  rb_define_module_function(mSpice, "furnsh", furnsh, -1);
  rb_define_module_function(mSpice, "gdpool", gdpool, 2);
  rb_define_module_function(mSpice, "getfov", getfov, 1);
  rb_define_module_function(mSpice, "gfdist", gfdist, -1);
  rb_define_module_function(mSpice, "gfsntc", gfsntc, -1);
  rb_define_module_function(mSpice, "gfsep", gfsep, -1);
  rb_define_module_function(mSpice, "gftfov", gftfov, -1);
  rb_define_module_function(mSpice, "gfrfov", gfrfov, -1);
  rb_define_module_function(mSpice, "ilumin", ilumin, -1);
  rb_define_module_function(mSpice, "kclear", kclear, 0);
  rb_define_module_function(mSpice, "ktotal", ktotal, -1);
  rb_define_module_function(mSpice, "latrec", latrec, 3);
  rb_define_module_function(mSpice, "lspcn", lspcn, 3);
  rb_define_module_function(mSpice, "m2eul", m2eul, 4);
  rb_define_module_function(mSpice, "m2q", m2q, -1);
  rb_define_module_function(mSpice, "mxv", mxv, 2);
  rb_define_module_function(mSpice, "gfoclt", gfoclt, -1);
  rb_define_module_function(mSpice, "pxform", pxform, 3);
  rb_define_module_function(mSpice, "q2m", q2m, -1);
  rb_define_module_function(mSpice, "sxform", sxform, 3);
  rb_define_module_function(mSpice, "reclat", reclat, -1);
  rb_define_module_function(mSpice, "recrad", recrad, -1);
  rb_define_module_function(mSpice, "rpd", rpd, 0);
  rb_define_module_function(mSpice, "sce2c", sce2c, 2);
  rb_define_module_function(mSpice, "scs2e", scs2e, 2);
  rb_define_module_function(mSpice, "sctiks", sctiks, 2);
  rb_define_module_function(mSpice, "sincpt", sincpt, -1);
  rb_define_module_function(mSpice, "spkcov", spkcov, -1);
  rb_define_module_function(mSpice, "spkezp", spkezp, -1);
  rb_define_module_function(mSpice, "spkezr", spkezr, -1);
  rb_define_module_function(mSpice, "spkpos", spkpos, 5);
  rb_define_module_function(mSpice, "str2et", str2et, -1);
  rb_define_module_function(mSpice, "subpnt", subpnt, -1);
  rb_define_module_function(mSpice, "subslr", subslr, -1);
  rb_define_module_function(mSpice, "ucrss", ucrss, 2);
  rb_define_module_function(mSpice, "unload", unload, -1);
  rb_define_module_function(mSpice, "vcrss", vcrss, 2);
  rb_define_module_function(mSpice, "vdist", vdist, 2);
  rb_define_module_function(mSpice, "vdot", vdot, 2);
  rb_define_module_function(mSpice, "vnorm", vnorm, -1);
  rb_define_module_function(mSpice, "vperp", vperp, 2);
  rb_define_module_function(mSpice, "vsep", vsep, 2);
  rb_define_module_function(mSpice, "vsub", vsub, 2);
#ifndef __APPLE__
  rb_define_module_function(mSpice, "wninsd",wninsd , -1);
#endif /*__APPLE__*/
  rb_define_module_function(mSpice, "xf2eul", xf2eul, 4);

  /*special ruby-only functions*/
  rb_define_module_function(mSpice, "et2rb", et2rb, 1);
  rb_define_module_function(mSpice, "rb2et", rb2et, 1);
  rb_define_module_function(mSpice, "wn2rb", wn2rb, 1);

  rb_eSpiceError = rb_define_class("SpiceError", rb_eStandardError);

  erract_c("SET", 0, (SpiceChar*)"RETURN");
  errprt_c("SET", 0, (SpiceChar*)"NONE");
}
