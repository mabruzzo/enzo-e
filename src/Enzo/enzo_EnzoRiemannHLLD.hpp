// Adapted from Athena++.
// See LICENSE_ATHENA++ file for license and copyright information

/// @file     enzo_EnzoRiemannHLLD.hpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     Thurs May 2 2019
/// @brief    [\ref Enzo] Enzo's HLLD approximate Riemann Solver.

// Currently, eint fluxes are computed by assuming that specific internal
// energy is a passive scalar. It may be worth considering the calculation of
// fluxes as though it's an actively advected quantity.

#ifndef ENZO_ENZO_RIEMANN_HLLD_HPP
#define ENZO_ENZO_RIEMANN_HLLD_HPP

typedef struct Cons1D {
  enzo_float d,mx,my,mz,e,by,bz;
} Cons1D;

#define SMALL_NUMBER 1.0e-8
#define SQR(x) (x*x)


struct HLLDImpl
{
  /// @class    EnzoRiemannHLLD
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Encapsulates operations of HLLD approximate Riemann
  /// Solver
public:
  static int scratch_space_length(const int n_cons_keys){
    return 2*n_cons_keys;
  }

  static void calc_riemann_fluxes
  (const enzo_float flux_l[], const enzo_float flux_r[],
   const enzo_float prim_l[], const enzo_float prim_r[],
   const enzo_float cons_l[], const enzo_float cons_r[],
   const enzo_float pressure_l, const enzo_float pressure_r,
   const EnzoAdvectionFieldLUT lut, const int n_keys,
   const bool barotropic_eos, const enzo_float gamma,
   const enzo_float isothermal_cs, const bool dual_energy,
   const int iz, const int iy, const int ix, EFlt3DArray flux_arrays[],
   enzo_float scratch_space[], enzo_float &vi_bar) throw()
  {

    enzo_float flxi[n_keys];             // temporary variable to store flux
    enzo_float spd[5];                   // signal speeds, left to right

    enzo_float gm1 = gamma - 1.0;

    Cons1D ul,ur;                  // L/R states, conserved variables (computed)
    Cons1D ulst,uldst,urdst,urst;  // Conserved variable for all states
    Cons1D fl,fr;                  // Fluxes for left & right states

    //--- Step 1.  Load L/R states into local variables

    // Preparation of prim_l and prim_r was handled ahead of time!
    const enzo_float* wli = prim_l;
    const enzo_float* wri = prim_r;

    enzo_float bxi = wli[lut.bfield_i];

    // Compute L/R states for selected conserved variables
    enzo_float bxsq = bxi*bxi;
    // (KGF): group transverse vector components for floating point
    // associativity symmetry
    enzo_float pbl = 0.5*(bxsq + (SQR(wli[lut.bfield_j]) +
				  SQR(wli[lut.bfield_k])));  
    enzo_float pbr = 0.5*(bxsq + (SQR(wri[lut.bfield_j]) +
				  SQR(wri[lut.bfield_k])));
    // kinetic energy (l/r)
    enzo_float kel = 0.5*wli[lut.density]*(SQR(wli[lut.velocity_i]) +
				   (SQR(wli[lut.velocity_j]) +
				    SQR(wli[lut.velocity_k])));
    enzo_float ker = 0.5*wri[lut.density]*(SQR(wri[lut.velocity_i]) +
				   (SQR(wri[lut.velocity_j]) +
				    SQR(wri[lut.velocity_k])));

    ul.d  = wli[lut.density];
    ul.mx = wli[lut.velocity_i]*ul.d;
    ul.my = wli[lut.velocity_j]*ul.d;
    ul.mz = wli[lut.velocity_k]*ul.d;
    ul.e  = pressure_l/gm1 + kel + pbl;
    ul.by = wli[lut.bfield_j];
    ul.bz = wli[lut.bfield_k];

    ur.d  = wri[lut.density];
    ur.mx = wri[lut.velocity_i]*ur.d;
    ur.my = wri[lut.velocity_j]*ur.d;
    ur.mz = wri[lut.velocity_k]*ur.d;
    ur.e  = pressure_r/gm1 + ker + pbr;
    ur.by = wri[lut.bfield_j];
    ur.bz = wri[lut.bfield_k];

    //--- Step 2.  Compute left & right wave speeds according to Miyoshi &
    // Kusano, eqn. (67)

    //enzo_float cfl = FastMagnetosonicSpeed(wli, bxi, pressure_l, lut, gamma);
    //enzo_float cfr = FastMagnetosonicSpeed(wri, bxi, pressure_r, lut, gamma);

    enzo_float cfl = EnzoRiemann::fast_magnetosonic_speed_(prim_l, lut,
							   pressure_l, gamma);
    enzo_float cfr = EnzoRiemann::fast_magnetosonic_speed_(prim_r, lut,
							   pressure_r, gamma);

    // This was originally included in the code
    //spd[0] = std::min( wli[lut.velocity_i]-cfl,
    //		       wri[lut.velocity_i]-cfr );
    //spd[4] = std::max( wli[lut.velocity_i]+cfl,
    //		       wri[lut.velocity_i]+cfr );

    
    // The following was originally commented out, but it is a more faithful
    // adaptation of Miyoshi & Kusano, eqn. (67)
    enzo_float cfmax = std::max(cfl,cfr);
    if (wli[lut.velocity_i] <= wri[lut.velocity_i]) {
      spd[0] = wli[lut.velocity_i] - cfmax;
      spd[4] = wri[lut.velocity_i] + cfmax;
    } else {
      spd[0] = wri[lut.velocity_i] - cfmax;
      spd[4] = wli[lut.velocity_i] + cfmax;
    }

    //--- Step 3.  Compute L/R fluxes

    enzo_float ptl = pressure_l + pbl; // total pressures L,R
    enzo_float ptr = pressure_r + pbr;

    fl.d  = ul.mx;
    fl.mx = ul.mx*wli[lut.velocity_i] + ptl - bxsq;
    fl.my = ul.my*wli[lut.velocity_i] - bxi*ul.by;
    fl.mz = ul.mz*wli[lut.velocity_i] - bxi*ul.bz;
    fl.e  = wli[lut.velocity_i]*(ul.e + ptl - bxsq) -
      bxi*(wli[lut.velocity_j]*ul.by + wli[lut.velocity_k]*ul.bz);
    fl.by = ul.by*wli[lut.velocity_i] - bxi*wli[lut.velocity_j];
    fl.bz = ul.bz*wli[lut.velocity_i] - bxi*wli[lut.velocity_k];

    fr.d  = ur.mx;
    fr.mx = ur.mx*wri[lut.velocity_i] + ptr - bxsq;
    fr.my = ur.my*wri[lut.velocity_i] - bxi*ur.by;
    fr.mz = ur.mz*wri[lut.velocity_i] - bxi*ur.bz;
    fr.e  = wri[lut.velocity_i]*(ur.e + ptr - bxsq) -
      bxi*(wri[lut.velocity_j]*ur.by + wri[lut.velocity_k]*ur.bz);
    fr.by = ur.by*wri[lut.velocity_i] - bxi*wri[lut.velocity_j];
    fr.bz = ur.bz*wri[lut.velocity_i] - bxi*wri[lut.velocity_k];

    //--- Step 4.  Compute middle and Alfven wave speeds

    enzo_float sdl = spd[0] - wli[lut.velocity_i];  // S_i-u_i (i=L or R)
    enzo_float sdr = spd[4] - wri[lut.velocity_i];

    // S_M: eqn (38) of Miyoshi & Kusano
    // (KGF): group ptl, ptr terms for floating point associativity symmetry
    spd[2] = (sdr*ur.mx - sdl*ul.mx + (ptl - ptr))/(sdr*ur.d - sdl*ul.d);

    enzo_float sdml   = spd[0] - spd[2];  // S_i-S_M (i=L or R)
    enzo_float sdmr   = spd[4] - spd[2];
    enzo_float sdml_inv = 1.0/sdml;
    enzo_float sdmr_inv = 1.0/sdmr;
    // eqn (43) of Miyoshi & Kusano
    ulst.d = ul.d * sdl * sdml_inv;
    urst.d = ur.d * sdr * sdmr_inv;
    enzo_float ulst_d_inv = 1.0/ulst.d;
    enzo_float urst_d_inv = 1.0/urst.d;
    enzo_float sqrtdl = std::sqrt(ulst.d);
    enzo_float sqrtdr = std::sqrt(urst.d);

    // eqn (51) of Miyoshi & Kusano
    spd[1] = spd[2] - fabs(bxi)/sqrtdl;
    spd[3] = spd[2] + fabs(bxi)/sqrtdr;

    //--- Step 5.  Compute intermediate states
    // eqn (23) explicitly becomes eq (41) of Miyoshi & Kusano
    // TODO(kfelker): place an assertion that ptstl==ptstr
    enzo_float ptstl = ptl + ul.d*sdl*(spd[2]-wli[lut.velocity_i]);
    enzo_float ptstr = ptr + ur.d*sdr*(spd[2]-wri[lut.velocity_i]);
    // these equations had issues when averaged
    // enzo_float ptstl = ptl + ul.d*sdl*(sdl-sdml);
    // enzo_float ptstr = ptr + ur.d*sdr*(sdr-sdmr);
    enzo_float ptst = 0.5*(ptstr + ptstl);  // total pressure (star state)

    // ul* - eqn (39) of M&K
    ulst.mx = ulst.d * spd[2];
    if (fabs(ul.d*sdl*sdml-bxsq) < (SMALL_NUMBER)*ptst) {
      // Degenerate case
      ulst.my = ulst.d * wli[lut.velocity_j];
      ulst.mz = ulst.d * wli[lut.velocity_k];

      ulst.by = ul.by;
      ulst.bz = ul.bz;
    } else {
      // eqns (44) and (46) of M&K
      enzo_float tmp = bxi*(sdl - sdml)/(ul.d*sdl*sdml - bxsq);
      ulst.my = ulst.d * (wli[lut.velocity_j] - ul.by*tmp);
      ulst.mz = ulst.d * (wli[lut.velocity_k] - ul.bz*tmp);

      // eqns (45) and (47) of M&K
      tmp = (ul.d*SQR(sdl) - bxsq)/(ul.d*sdl*sdml - bxsq);
      ulst.by = ul.by * tmp;
      ulst.bz = ul.bz * tmp;
    }
    // v_i* dot B_i*
    // (KGF): group transverse momenta terms for floating point
    // associativity symmetry
    enzo_float vbstl = (ulst.mx*bxi+(ulst.my*ulst.by+ulst.mz*ulst.bz))*ulst_d_inv;
    // eqn (48) of M&K
    // (KGF): group transverse by, bz terms for floating point
    // associativity symmetry
    ulst.e = (sdl*ul.e - ptl*wli[lut.velocity_i] + ptst*spd[2] +
              bxi*(wli[lut.velocity_i]*bxi + (wli[lut.velocity_j]*ul.by + wli[lut.velocity_k]*ul.bz) - vbstl))*sdml_inv;

  // ur* - eqn (39) of M&K
    urst.mx = urst.d * spd[2];
    if (fabs(ur.d*sdr*sdmr - bxsq) < (SMALL_NUMBER)*ptst) {
      // Degenerate case
      urst.my = urst.d * wri[lut.velocity_j];
      urst.mz = urst.d * wri[lut.velocity_k];

      urst.by = ur.by;
      urst.bz = ur.bz;
    } else {
      // eqns (44) and (46) of M&K
      enzo_float tmp = bxi*(sdr - sdmr)/(ur.d*sdr*sdmr - bxsq);
      urst.my = urst.d * (wri[lut.velocity_j] - ur.by*tmp);
      urst.mz = urst.d * (wri[lut.velocity_k] - ur.bz*tmp);

      // eqns (45) and (47) of M&K
      tmp = (ur.d*SQR(sdr) - bxsq)/(ur.d*sdr*sdmr - bxsq);
      urst.by = ur.by * tmp;
      urst.bz = ur.bz * tmp;
    }
    // v_i* dot B_i*
    // (KGF): group transverse momenta terms for floating point associativity symmetry
    enzo_float vbstr = (urst.mx*bxi+(urst.my*urst.by+urst.mz*urst.bz))*urst_d_inv;
    // eqn (48) of M&K
    // (KGF): group transverse by, bz terms for floating point associativity symmetry
    urst.e = (sdr*ur.e - ptr*wri[lut.velocity_i] + ptst*spd[2] +
              bxi*(wri[lut.velocity_i]*bxi + (wri[lut.velocity_j]*ur.by + wri[lut.velocity_k]*ur.bz) - vbstr))*sdmr_inv;
    // ul** and ur** - if Bx is near zero, same as *-states
    if (0.5*bxsq < (SMALL_NUMBER)*ptst) {
      uldst = ulst;
      urdst = urst;
    } else {
      enzo_float invsumd = 1.0/(sqrtdl + sqrtdr);
      enzo_float bxsig = (bxi > 0.0 ? 1.0 : -1.0);

      uldst.d = ulst.d;
      urdst.d = urst.d;

      uldst.mx = ulst.mx;
      urdst.mx = urst.mx;

      // eqn (59) of M&K
      enzo_float tmp = invsumd*(sqrtdl*(ulst.my*ulst_d_inv) + sqrtdr*(urst.my*urst_d_inv) +
                 bxsig*(urst.by - ulst.by));
      uldst.my = uldst.d * tmp;
      urdst.my = urdst.d * tmp;

      // eqn (60) of M&K
      tmp = invsumd*(sqrtdl*(ulst.mz*ulst_d_inv) + sqrtdr*(urst.mz*urst_d_inv) +
            bxsig*(urst.bz - ulst.bz));
      uldst.mz = uldst.d * tmp;
      urdst.mz = urdst.d * tmp;

      // eqn (61) of M&K
      tmp = invsumd*(sqrtdl*urst.by + sqrtdr*ulst.by +
                     bxsig*sqrtdl*sqrtdr*((urst.my*urst_d_inv) - (ulst.my*ulst_d_inv)));
      uldst.by = urdst.by = tmp;

      // eqn (62) of M&K
      tmp = invsumd*(sqrtdl*urst.bz + sqrtdr*ulst.bz +
                     bxsig*sqrtdl*sqrtdr*((urst.mz*urst_d_inv) - (ulst.mz*ulst_d_inv)));
      uldst.bz = urdst.bz = tmp;

      // eqn (63) of M&K
      tmp = spd[2]*bxi + (uldst.my*uldst.by + uldst.mz*uldst.bz)/uldst.d;
      uldst.e = ulst.e - sqrtdl*bxsig*(vbstl - tmp);
      urdst.e = urst.e + sqrtdr*bxsig*(vbstr - tmp);
    }

//--- Step 6.  Compute flux
    uldst.d = spd[1] * (uldst.d - ulst.d);
    uldst.mx = spd[1] * (uldst.mx - ulst.mx);
    uldst.my = spd[1] * (uldst.my - ulst.my);
    uldst.mz = spd[1] * (uldst.mz - ulst.mz);
    uldst.e = spd[1] * (uldst.e - ulst.e);
    uldst.by = spd[1] * (uldst.by - ulst.by);
    uldst.bz = spd[1] * (uldst.bz - ulst.bz);

    ulst.d = spd[0] * (ulst.d - ul.d);
    ulst.mx = spd[0] * (ulst.mx - ul.mx);
    ulst.my = spd[0] * (ulst.my - ul.my);
    ulst.mz = spd[0] * (ulst.mz - ul.mz);
    ulst.e = spd[0] * (ulst.e - ul.e);
    ulst.by = spd[0] * (ulst.by - ul.by);
    ulst.bz = spd[0] * (ulst.bz - ul.bz);

    urdst.d = spd[3] * (urdst.d - urst.d);
    urdst.mx = spd[3] * (urdst.mx - urst.mx);
    urdst.my = spd[3] * (urdst.my - urst.my);
    urdst.mz = spd[3] * (urdst.mz - urst.mz);
    urdst.e = spd[3] * (urdst.e - urst.e);
    urdst.by = spd[3] * (urdst.by - urst.by);
    urdst.bz = spd[3] * (urdst.bz - urst.bz);

    urst.d = spd[4] * (urst.d  - ur.d);
    urst.mx = spd[4] * (urst.mx - ur.mx);
    urst.my = spd[4] * (urst.my - ur.my);
    urst.mz = spd[4] * (urst.mz - ur.mz);
    urst.e = spd[4] * (urst.e - ur.e);
    urst.by = spd[4] * (urst.by - ur.by);
    urst.bz = spd[4] * (urst.bz - ur.bz);

    if (spd[0] >= 0.0) {
      // return Fl if flow is supersonic
      flxi[lut.density] = fl.d;
      flxi[lut.velocity_i] = fl.mx;
      flxi[lut.velocity_j] = fl.my;
      flxi[lut.velocity_k] = fl.mz;
      flxi[lut.total_energy] = fl.e;
      flxi[lut.bfield_j] = fl.by;
      flxi[lut.bfield_k] = fl.bz;

      vi_bar = prim_l[lut.velocity_i];
    } else if (spd[4] <= 0.0) {
      // return Fr if flow is supersonic
      flxi[lut.density] = fr.d;
      flxi[lut.velocity_i] = fr.mx;
      flxi[lut.velocity_j] = fr.my;
      flxi[lut.velocity_k] = fr.mz;
      flxi[lut.total_energy] = fr.e;
      flxi[lut.bfield_j] = fr.by;
      flxi[lut.bfield_k] = fr.bz;

      vi_bar = prim_r[lut.velocity_i];
    } else if (spd[1] >= 0.0) {
      // return Fl*
      flxi[lut.density] = fl.d  + ulst.d;
      flxi[lut.velocity_i] = fl.mx + ulst.mx;
      flxi[lut.velocity_j] = fl.my + ulst.my;
      flxi[lut.velocity_k] = fl.mz + ulst.mz;
      flxi[lut.total_energy] = fl.e  + ulst.e;
      flxi[lut.bfield_j] = fl.by + ulst.by;
      flxi[lut.bfield_k] = fl.bz + ulst.bz;

      vi_bar = spd[2] * sdl/sdml;
    } else if (spd[2] >= 0.0) {
      // return Fl**
      flxi[lut.density] = fl.d  + ulst.d + uldst.d;
      flxi[lut.velocity_i] = fl.mx + ulst.mx + uldst.mx;
      flxi[lut.velocity_j] = fl.my + ulst.my + uldst.my;
      flxi[lut.velocity_k] = fl.mz + ulst.mz + uldst.mz;
      flxi[lut.total_energy] = fl.e  + ulst.e + uldst.e;
      flxi[lut.bfield_j] = fl.by + ulst.by + uldst.by;
      flxi[lut.bfield_k] = fl.bz + ulst.bz + uldst.bz;

      vi_bar = spd[2] * sdl/sdml;
    } else if (spd[3] > 0.0) {
      // return Fr**
      flxi[lut.density] = fr.d + urst.d + urdst.d;
      flxi[lut.velocity_i] = fr.mx + urst.mx + urdst.mx;
      flxi[lut.velocity_j] = fr.my + urst.my + urdst.my;
      flxi[lut.velocity_k] = fr.mz + urst.mz + urdst.mz;
      flxi[lut.total_energy] = fr.e + urst.e + urdst.e;
      flxi[lut.bfield_j] = fr.by + urst.by + urdst.by;
      flxi[lut.bfield_k] = fr.bz + urst.bz + urdst.bz;

      vi_bar = spd[2] * sdr/sdmr;
    } else {
      // return Fr*
      flxi[lut.density] = fr.d  + urst.d;
      flxi[lut.velocity_i] = fr.mx + urst.mx;
      flxi[lut.velocity_j] = fr.my + urst.my;
      flxi[lut.velocity_k] = fr.mz + urst.mz;
      flxi[lut.total_energy] = fr.e  + urst.e;
      flxi[lut.bfield_j] = fr.by + urst.by;
      flxi[lut.bfield_k] = fr.bz + urst.bz;

      vi_bar = spd[2] * sdr/sdmr;
    }

    if (dual_energy) {
      if (spd[0] >= 0.0 || spd[1] >= 0.0 || spd[2] >= 0.0) {
	flux_arrays[lut.internal_energy](iz,iy,ix)
	  = prim_l[lut.internal_energy] * flxi[lut.density];
      } else {
	flux_arrays[lut.internal_energy](iz,iy,ix)
	  = prim_r[lut.internal_energy] * flxi[lut.density];
      }
    }

    flux_arrays[lut.density](iz,iy,ix) = flxi[lut.density];
    flux_arrays[lut.velocity_i](iz,iy,ix) = flxi[lut.velocity_i];
    flux_arrays[lut.velocity_j](iz,iy,ix) = flxi[lut.velocity_j];
    flux_arrays[lut.velocity_k](iz,iy,ix) = flxi[lut.velocity_k];
    flux_arrays[lut.total_energy](iz,iy,ix) = flxi[lut.total_energy];
    // The following is handled slightly differently from Athena++
    flux_arrays[lut.bfield_i](iz,iy,ix) = 0.0;
    flux_arrays[lut.bfield_j](iz,iy,ix) = flxi[lut.bfield_j];
    flux_arrays[lut.bfield_k](iz,iy,ix) = flxi[lut.bfield_k];
  }
};


/// @class    EnzoRiemannHLLD
/// @ingroup  Enzo
/// @brief    [\ref Enzo] Encapsulates HLLD approximate Riemann Solver
using EnzoRiemannHLLD = EnzoRiemannImpl<HLLDImpl>;


#endif /* ENZO_ENZO_RIEMANN_HLLD_HPP */
