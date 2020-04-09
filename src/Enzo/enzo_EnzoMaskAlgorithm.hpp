// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMashAlgorithm.hpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     Tues March 10 2020
/// @brief    [\ref Enzo] Definition of algorithms for iterating over elements
/// in boolean CellArray instances. This is mainly employed to implement
/// Riemann Solver Fallback

template<class InnerLoopFunc>
inline void mask_iter_helper_(CelloArray<bool,3> &mask,
			      InnerLoopFunc inner_loop_func,
			      CelloArray<bool,3> &visited, int step,
			      int delta_ind[3])
{
  int max_upper[3]; // important to set this with visited and not mask
  for (int i = 0; i < 3; i++){ max_upper[i] = visited.shape(i); }

  // NOTE: mask can have a different shape from visited
  for (int iz=0; iz<mask.shape(0); iz++) {
    for (int iy=0; iy<mask.shape(1); iy++) {
      for (int ix=0; ix<mask.shape(2); ix++) {
	if (mask(iz,iy,ix)){
	  int cur[3];    cur[0] = iz;    cur[1] = iy;    cur[2] = ix;
	  int lower[3]; // kl, jl, il
	  int upper[3]; // ku, ju, iu
	  for (int i=0; i < 3; i++){
	    lower[i] = std::max(cur[i] - delta_ind[i],0);
	    upper[i] = std::max(cur[i] + 1 + step*delta_ind[i], max_upper[i]);
	  }

	  // adjust lower iteration limits to account for overlap with past
	  // iteration and update visited so that we don't overlap in future
	  for (int iz2 = lower[0]; iz2 < upper[0]; iz2++){
	    for (int iy2 = lower[1]; iy2 < upper[1]; iy2++){
	      for (int ix2 = lower[2]; ix2 < upper[2]; ix2++){
		if (!visited(iz2,iy2,ix2)){ break; }
		for (int i = 0; i<3; i++){lower[i] += delta_ind[i];}
		visited(iz2,iy2,ix2) = true;
	      }
	    }
	  }
	  inner_loop_func(lower[0], upper[0], lower[1], upper[1],
			  lower[2], upper[2]);
	}
      }
    }
  }

}

/// Metafunction that applies a functor at every cell-centered value specified
/// by a boolean mask
///
/// @details
/// When `line_length` is `1`, then this literally just executes
/// `inner_loop_func` at every location, `(iz,iy,ix)`, where `mask(iz,iy,ix)`.
/// Otherwise, the `inner_loop_func` is evaluated at all of these locations
/// plus the `line_length//2` neigboring indices along dimension `dim` (note
/// that `line_length` must be odd)
///
/// An equivalent description is that this evaluates the functor at all
/// locations where the convolution of `mask` and a kernel, specified by the
/// arguments, evaluates to `true`. The kernel is effectively a 3D array of
/// true values with `line_length` elements along dimension `dim` and just 1
/// element along the other two axes - it is effectively a line. The kernel is
/// centered on the element at `line_length//2 + 1`
///
/// @tparam InnerLoopFunc type of the functor to call.
/// @param mask (partially) specifies locations to be iterated over
/// @param inner_loop_func functor to evaluate at each location. This must have
///     the function signature
///     `void fun(int kl, int ku, int jl, int ju, int il, int iu);` where
///     `kl,ku`, `jl,ju`, and `il,iu` represent pairs of starting and stopping
///     indices to iterate over along the z, y, and x axes.
/// @param dim the functor will be evaluated at neigbors `true` locations in
///     the mask along this dimension.
/// @param line_length specifies the length of the convolution kernel along
///     dimension `dim` that is effectively used to determine the locations
///     where inner_loop_functor gets evaluated. This must be odd and positive.
///
/// @par Note
/// This can be optimized to avoid allocating an array to track locations that
/// have been visited for the cases when `line_length` is `3` (arguably the 
/// most commonly used case) or `1`
template<class InnerLoopFunc>
void mask_iter(CelloArray<bool,3> &mask, InnerLoopFunc inner_loop_func,
	       int dim, int line_length)
{  
  ASSERT("mask_iter", "line_length must be odd and positive",
	 line_length > 3 && (line_length % 2 == 1));
  // determine how many elements should be visited to the left of line elements
  int line_neighbors = (line_length - 1) / 2;

  EnzoPermutedCoordinates coord(dim);
  int delta_ind[3]; // = {delta_iz, delta_iy, delta_ix}
  coord.i_unit_vector(delta_ind[2], delta_ind[1], delta_ind[0]);

  // create an array to keep track of locations that have already been visited
  // all values are initialized to false, by default
  CelloArray<bool,3> visited(mask.shape(0), mask.shape(1), mask.shape(2));

  mask_iter_helper_(mask, inner_loop_func, visited, line_neighbors, delta_ind);
}

/// Evaluates a functor at all cell-interface locations (centered along the
/// `dim`-cell faces) immediately adjacent to the (cell-centered) locations in
/// mask that are set to `true`
///
/// @details
/// This function assumes that the allowed locations where `inner_loop_func`
/// can be evaluated include the faces exterior to `mask`. Relatedly, it
/// assumes that a value, `index`, it passes to `inner_loop_func` as an
/// iteration limit along dimension `dim` corresponds to the cell-interface
/// at `index-0.5`.
///
/// To be clear about the consequences of the above assumptions, consider the
/// bounds that would be passed to a functor that would achieve an equivalent
/// result to using this function to evaluate the functor with a `mask` that is
/// `true` everywhere. If the `mask` has `n` elements along axis `dim`, then
/// the lower and upper bounds for this dimension that need to be passed to the
/// functor are `0` and `n+1`. For the other dimensions the lower bounds should
/// be `0` while the upper bounds should be the number of elements along the
/// corresponding axis of `mask`.
///
/// @tparam InnerLoopFunc type of the functor to call.
/// @param mask specifies the cell-centered masked locations.
/// @param inner_loop_func functor to evaluate at each location. This must have
///     the function signature
///     `void fun(int kl, int ku, int jl, int ju, int il, int iu);` where
///     `kl,ku`, `jl,ju`, and `il,iu` represent pairs of starting and stopping
///     indices to iterate over along the z, y, and x axes.
/// @param dim Indicates the dimension along which the functor treats the
///     iteration indices as face-centered (e.g. if dim == 0, then the
///     dimension, then a given index `i` actually corresponds to the face at
///     `i-1/2`).
///
/// @par Note
/// This can be optimized to avoid allocating an array to track locations that
/// have been visited.
template<class InnerLoopFunc>
void mask_neighbor_face_iter(CelloArray<bool,3> &mask,
			     InnerLoopFunc inner_loop_func, int dim){
  EnzoPermutedCoordinates coord(dim);
  int delta_ind[3]; // = {delta_iz, delta_iy, delta_ix}
  coord.i_unit_vector(delta_ind[2], delta_ind[1], delta_ind[0]);
  // create an array to keep track of locations that have already been visited
  // all values are initialized to false, by default
  CelloArray<bool,3> visited(mask.shape(0) + delta_ind[0],
			     mask.shape(1) + delta_ind[1],
			     mask.shape(2) + delta_ind[2]);
  mask_iter_helper_(mask, inner_loop_func, visited, 1, delta_ind);
}
