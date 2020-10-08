#include "cello.hpp"
#include "enzo.hpp"

//----------------------------------------------------------------------

const EFlt3DArray& EnzoEFltArrayMap::at(const std::string& key) const noexcept
{
  auto result = map_.find(key);
  if (result == map_.cend()){
    ERROR1("EnzoEFltArrayMap::at", "map doesn't contain the key: \"%s\"",
           key.c_str());
  } else {
    return result->second;
  }
}

//----------------------------------------------------------------------

const EFlt3DArray exclude_stale_cells_(const EFlt3DArray &arr, int stale_depth)
{
  ASSERT("exclude_stale_cells_","each dim of arr must exceed 2*stale_depth.",
	 2*stale_depth < arr.shape(0) && 2*stale_depth < arr.shape(1) &&
	 2*stale_depth < arr.shape(2));

  return arr.subarray(CSlice(stale_depth, -stale_depth),
		      CSlice(stale_depth, -stale_depth),
		      CSlice(stale_depth, -stale_depth));
}


const EFlt3DArray EnzoEFltArrayMap::get(const std::string& key,
                                        int stale_depth) const noexcept
{
  ASSERT("EnzoEFltArrayMap::get", "stale_depth must be >= 0",
         stale_depth >= 0);
  const EFlt3DArray& arr = this->at(key);
  if (stale_depth > 0){
    return exclude_stale_cells_(arr,stale_depth);
  } else {
    return arr;
  }
}

