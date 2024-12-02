#ifndef BLOCKINFO_HPP
#define BLOCKINFO_HPP

#include "Geometry.hpp"
#include <map>

struct BlockInfo_STRUCT {
  local_int_t start;
  local_int_t end;
};
typedef struct BlockInfo_STRUCT BlockInfo;

struct BlockLayout_STRUCT {
  local_int_t numBlocks;
  std::vector<BlockInfo> blkInfo;
  std::vector< std::vector<local_int_t> > blkDepList;
  std::vector< std::vector<local_int_t> > colorBlkInfo;
  std::vector< std::vector<local_int_t> > coarserToCurrentBlkInfo;
  // Each neighbor will have for each block a list of indices to send
  local_int_t extraBlocks;
  std::vector< std::map<local_int_t, std::vector<local_int_t> > > neighborHaloToBlkHalo;
};
typedef struct BlockLayout_STRUCT BlockLayout;

#endif // BLOCKINFO_HPP

