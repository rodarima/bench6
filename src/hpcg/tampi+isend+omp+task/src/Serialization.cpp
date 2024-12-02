#include <fstream>
#include <cstddef>

#include "SparseMatrix.hpp"
#include "Geometry.hpp"
#include "hpcg.hpp"
#include "Serialization.hpp"

void ParamsSerialize(HPCG_Params &params, const std::string &fname) {
  std::ofstream ofile;
  ofile.open(fname.c_str(), std::ios::out | std::ios::binary | std::ios::trunc);
  die_if(ofile.fail(), "HPCG_Params: open() error");
  ofile.write(reinterpret_cast<char*>(&params), sizeof(params));
  ofile.close();
}

void ParamsDeserialize(HPCG_Params &params, const std::string &fname) {
  std::ifstream ifile;
  ifile.open(fname.c_str(), std::ios::in | std::ios::binary);
  die_if(ifile.fail(), "HPCG_Params: open() error");
  ifile.read(reinterpret_cast<char*>(&params), sizeof(params));
  // We do not care about load/store paths when deserializing.
  // Set pointers to nullptr to avoid freeng on DeleteParamsData.
  params.load_path = 0;
  params.store_path = 0;
  ifile.close();
}

void GeometrySerialize(Geometry &geom, const std::string &fname) {
  std::ofstream ofile;
  ofile.open(fname.c_str(), std::ios::out | std::ios::binary | std::ios::trunc);
  die_if(ofile.fail(), "Geometry: open() error");
  ofile.write(reinterpret_cast<char*>(&geom.size), sizeof(geom.size));
  ofile.write(reinterpret_cast<char*>(&geom.rank), sizeof(geom.rank));
  ofile.write(reinterpret_cast<char*>(&geom.numThreads), sizeof(geom.numThreads));
  ofile.write(reinterpret_cast<char*>(&geom.nx), sizeof(geom.nx));
  ofile.write(reinterpret_cast<char*>(&geom.ny), sizeof(geom.ny));
  ofile.write(reinterpret_cast<char*>(&geom.nz), sizeof(geom.nz));
  ofile.write(reinterpret_cast<char*>(&geom.npx), sizeof(geom.npx));
  ofile.write(reinterpret_cast<char*>(&geom.npy), sizeof(geom.npy));
  ofile.write(reinterpret_cast<char*>(&geom.npz), sizeof(geom.npz));
  ofile.write(reinterpret_cast<char*>(&geom.pz), sizeof(geom.pz));
  ofile.write(reinterpret_cast<char*>(&geom.npartz), sizeof(geom.npartz));
  if (geom.pz == 0) {
    ofile.write(reinterpret_cast<char*>(&geom.partz_ids[0]), sizeof(int[1]));
    ofile.write(reinterpret_cast<char*>(&geom.partz_nz[0]), sizeof(local_int_t[1]));
  }
  else {
    ofile.write(reinterpret_cast<char*>(&geom.partz_ids[0]), sizeof(int[2]));
    ofile.write(reinterpret_cast<char*>(&geom.partz_nz[0]), sizeof(local_int_t[2]));
  }
  ofile.write(reinterpret_cast<char*>(&geom.ipx), sizeof(geom.ipx));
  ofile.write(reinterpret_cast<char*>(&geom.ipy), sizeof(geom.ipy));
  ofile.write(reinterpret_cast<char*>(&geom.ipz), sizeof(geom.ipz));
  ofile.write(reinterpret_cast<char*>(&geom.gnx), sizeof(geom.gnx));
  ofile.write(reinterpret_cast<char*>(&geom.gny), sizeof(geom.gny));
  ofile.write(reinterpret_cast<char*>(&geom.gnz), sizeof(geom.gnz));
  ofile.write(reinterpret_cast<char*>(&geom.gix0), sizeof(geom.gix0));
  ofile.write(reinterpret_cast<char*>(&geom.giy0), sizeof(geom.giy0));
  ofile.write(reinterpret_cast<char*>(&geom.giz0), sizeof(geom.giz0));
  ofile.close();
}

void GeometryDeserialize(Geometry &geom, const std::string &fname) {
  std::ifstream ifile;
  ifile.open(fname.c_str(), std::ios::in | std::ios::binary);
  die_if(ifile.fail(), "Geometry: open() error");
  ifile.read(reinterpret_cast<char*>(&geom.size), sizeof(geom.size));
  ifile.read(reinterpret_cast<char*>(&geom.rank), sizeof(geom.rank));
  ifile.read(reinterpret_cast<char*>(&geom.numThreads), sizeof(geom.numThreads));
  ifile.read(reinterpret_cast<char*>(&geom.nx), sizeof(geom.nx));
  ifile.read(reinterpret_cast<char*>(&geom.ny), sizeof(geom.ny));
  ifile.read(reinterpret_cast<char*>(&geom.nz), sizeof(geom.nz));
  ifile.read(reinterpret_cast<char*>(&geom.npx), sizeof(geom.npx));
  ifile.read(reinterpret_cast<char*>(&geom.npy), sizeof(geom.npy));
  ifile.read(reinterpret_cast<char*>(&geom.npz), sizeof(geom.npz));
  ifile.read(reinterpret_cast<char*>(&geom.pz), sizeof(geom.pz));
  ifile.read(reinterpret_cast<char*>(&geom.npartz), sizeof(geom.npartz));
  if (geom.pz == 0) {
    geom.partz_ids = new int[1];
    geom.partz_nz = new local_int_t[1];
    ifile.read(reinterpret_cast<char*>(&geom.partz_ids[0]), sizeof(int[1]));
    ifile.read(reinterpret_cast<char*>(&geom.partz_nz[0]), sizeof(local_int_t[1]));
  }
  else {
    geom.partz_ids = new int[2];
    geom.partz_nz = new local_int_t[2];
    ifile.read(reinterpret_cast<char*>(&geom.partz_ids[0]), sizeof(int[2]));
    ifile.read(reinterpret_cast<char*>(&geom.partz_nz[0]), sizeof(local_int_t[2]));
  }
  ifile.read(reinterpret_cast<char*>(&geom.ipx), sizeof(geom.ipx));
  ifile.read(reinterpret_cast<char*>(&geom.ipy), sizeof(geom.ipy));
  ifile.read(reinterpret_cast<char*>(&geom.ipz), sizeof(geom.ipz));
  ifile.read(reinterpret_cast<char*>(&geom.gnx), sizeof(geom.gnx));
  ifile.read(reinterpret_cast<char*>(&geom.gny), sizeof(geom.gny));
  ifile.read(reinterpret_cast<char*>(&geom.gnz), sizeof(geom.gnz));
  ifile.read(reinterpret_cast<char*>(&geom.gix0), sizeof(geom.gix0));
  ifile.read(reinterpret_cast<char*>(&geom.giy0), sizeof(geom.giy0));
  ifile.read(reinterpret_cast<char*>(&geom.giz0), sizeof(geom.giz0));
  ifile.close();
}

void SparseMatrixSerialize(SparseMatrix &A, const std::string &fname) {
  std::ofstream ofile;
  ofile.open(fname.c_str(), std::ios::out | std::ios::binary | std::ios::trunc);
  die_if(ofile.fail(), "open() error");
  // write_params.write(reinterpret_cast<char*>(&params.title), sizeof(params.title));
  // write_params.write(reinterpret_cast<char*>(&params.geom), sizeof(params.geom));
  ofile.write(reinterpret_cast<char*>(&A.totalNumberOfRows), sizeof(A.totalNumberOfRows));
  ofile.write(reinterpret_cast<char*>(&A.totalNumberOfNonzeros), sizeof(A.totalNumberOfNonzeros));
  ofile.write(reinterpret_cast<char*>(&A.localNumberOfRows), sizeof(A.localNumberOfRows));
  ofile.write(reinterpret_cast<char*>(&A.localNumberOfColumns), sizeof(A.localNumberOfColumns));
  ofile.write(reinterpret_cast<char*>(&A.localNumberOfNonzeros), sizeof(A.localNumberOfNonzeros));
  ofile.write(&A.nonzerosInRow[0], A.localNumberOfRows*sizeof(A.nonzerosInRow[0]));
#ifndef HPCG_CONTIGUOUS_ARRAYS
  for (local_int_t i = 0; i < A.localNumberOfRows; i++) {
    int cur_nnz = A.nonzerosInRow[i];
    ofile.write(reinterpret_cast<char*>(&A.mtxIndG[i][0]), cur_nnz*sizeof(A.mtxIndG[0][0]));
  }
  for (local_int_t i = 0; i < A.localNumberOfRows; i++) {
    int cur_nnz = A.nonzerosInRow[i];
    ofile.write(reinterpret_cast<char*>(&A.mtxIndL[i][0]), cur_nnz*sizeof(A.mtxIndL[0][0]));
  }
  for (local_int_t i = 0; i < A.localNumberOfRows; i++) {
    int cur_nnz = A.nonzerosInRow[i];
    ofile.write(reinterpret_cast<char*>(&A.matrixValues[i][0]), cur_nnz*sizeof(A.matrixValues[0][0]));
  }
#else
  const local_int_t numberOfNonzerosPerRow = 27; // We are approximating a 27-point finite element/volume/offseterence 3D stencil
  ofile.write(reinterpret_cast<char*>(&A.mtxIndG[0][0]), A.localNumberOfRows*numberOfNonzerosPerRow*sizeof(A.mtxIndG[0][0]));
  ofile.write(reinterpret_cast<char*>(&A.mtxIndL[0][0]), A.localNumberOfRows*numberOfNonzerosPerRow*sizeof(A.mtxIndL[0][0]));
  ofile.write(reinterpret_cast<char*>(&A.matrixValues[0][0]), A.localNumberOfRows*numberOfNonzerosPerRow*sizeof(A.matrixValues[0][0]));
#endif
  // matrixDiagonal points to matrixValues, store offsets
  for (local_int_t i = 0; i < A.localNumberOfRows; i++) {
    ptrdiff_t offset = &A.matrixValues[i][0] - A.matrixDiagonal[i];
    ofile.write(reinterpret_cast<char*>(&offset), sizeof(offset));
  }
  // for (GlobalToLocalMap::iterator it = A.globalToLocalMap.begin();
  //      it != A.globalToLocalMap.end(); ++it) {
  //   ofile.write(reinterpret_cast<char*>(const_cast<global_int_t*>(&it->first)), sizeof(it->first));
  //   ofile.write(reinterpret_cast<char*>(&it->second), sizeof(it->second));
  // }
  ofile.write(reinterpret_cast<char*>(&A.localToGlobalMap[0]), A.localToGlobalMap.size()*sizeof(A.localToGlobalMap[0]));
  ofile.write(reinterpret_cast<char*>(&A.isDotProductOptimized), sizeof(A.isDotProductOptimized));
  ofile.write(reinterpret_cast<char*>(&A.isSpmvOptimized), sizeof(A.isSpmvOptimized));
  ofile.write(reinterpret_cast<char*>(&A.isMgOptimized), sizeof(A.isMgOptimized));
  ofile.write(reinterpret_cast<char*>(&A.isWaxpbyOptimized), sizeof(A.isWaxpbyOptimized));
  ofile.write(reinterpret_cast<char*>(&A.numberOfExternalValues), sizeof(A.numberOfExternalValues));
  ofile.write(reinterpret_cast<char*>(&A.numberOfSendNeighbors), sizeof(A.numberOfSendNeighbors));
  ofile.write(reinterpret_cast<char*>(&A.totalToBeSent), sizeof(A.totalToBeSent));
  ofile.write(reinterpret_cast<char*>(&A.elementsToSend[0]), A.totalToBeSent*sizeof(A.elementsToSend[0]));
  ofile.write(reinterpret_cast<char*>(&A.neighbors[0]), A.numberOfSendNeighbors*sizeof(A.neighbors[0]));
  ofile.write(reinterpret_cast<char*>(&A.receiveLength[0]), A.numberOfSendNeighbors*sizeof(A.receiveLength[0]));
  ofile.write(reinterpret_cast<char*>(&A.sendLength[0]), A.numberOfSendNeighbors*sizeof(A.sendLength[0]));
  // A.sendBuffer is just a tmp buffer
  ofile.close();
}

void SparseMatrixDeserialize(SparseMatrix &A, const std::string &fname) {
  std::ifstream ifile;
  ifile.open(fname.c_str(), std::ios::in | std::ios::binary);
  die_if(ifile.fail(), "SparseMatrix: open() error");
  // read_params.read(reinterpret_cast<char*>(&params.title), sizeof(params.title));
  // read_params.read(reinterpret_cast<char*>(&params.geom), sizeof(params.geom));
  ifile.read(reinterpret_cast<char*>(&A.totalNumberOfRows), sizeof(A.totalNumberOfRows));
  ifile.read(reinterpret_cast<char*>(&A.totalNumberOfNonzeros), sizeof(A.totalNumberOfNonzeros));
  ifile.read(reinterpret_cast<char*>(&A.localNumberOfRows), sizeof(A.localNumberOfRows));
  ifile.read(reinterpret_cast<char*>(&A.localNumberOfColumns), sizeof(A.localNumberOfColumns));
  ifile.read(reinterpret_cast<char*>(&A.localNumberOfNonzeros), sizeof(A.localNumberOfNonzeros));
  A.nonzerosInRow = new char[A.localNumberOfRows];
  A.mtxIndG = new global_int_t*[A.localNumberOfRows];
  A.mtxIndL = new local_int_t*[A.localNumberOfRows];
  A.matrixValues = new double*[A.localNumberOfRows];
  A.matrixDiagonal = new double*[A.localNumberOfRows];
  ifile.read(&A.nonzerosInRow[0], A.localNumberOfRows*sizeof(A.nonzerosInRow[0]));
#ifndef HPCG_CONTIGUOUS_ARRAYS
  const local_int_t numberOfNonzerosPerRow = 27; // We are approximating a 27-point finite element/volume/offseterence 3D stencil
  // Now allocate the arrays pointed to
  for (local_int_t i=0; i< A.localNumberOfRows; ++i)
    A.mtxIndL[i] = new local_int_t[numberOfNonzerosPerRow];
  for (local_int_t i=0; i< A.localNumberOfRows; ++i)
    A.matrixValues[i] = new double[numberOfNonzerosPerRow];
  for (local_int_t i=0; i< A.localNumberOfRows; ++i)
   A.mtxIndG[i] = new global_int_t[numberOfNonzerosPerRow];

  for (local_int_t i = 0; i < A.localNumberOfRows; i++) {
    int cur_nnz = A.nonzerosInRow[i];
    ifile.read(reinterpret_cast<char*>(&A.mtxIndG[i][0]), cur_nnz*sizeof(A.mtxIndG[0][0]));
  }
  for (local_int_t i = 0; i < A.localNumberOfRows; i++) {
    int cur_nnz = A.nonzerosInRow[i];
    ifile.read(reinterpret_cast<char*>(&A.mtxIndL[i][0]), cur_nnz*sizeof(A.mtxIndL[0][0]));
  }
  for (local_int_t i = 0; i < A.localNumberOfRows; i++) {
    int cur_nnz = A.nonzerosInRow[i];
    ifile.read(reinterpret_cast<char*>(&A.matrixValues[i][0]), cur_nnz*sizeof(A.matrixValues[0][0]));
  }

#else
  // Now allocate the arrays pointed to
  A.mtxIndL[0] = new local_int_t[A.localNumberOfRows * numberOfNonzerosPerRow];
  A.matrixValues[0] = new double[A.localNumberOfRows * numberOfNonzerosPerRow];
  A.mtxIndG[0] = new global_int_t[A.localNumberOfRows * numberOfNonzerosPerRow];

  for (local_int_t i=1; i< A.localNumberOfRows; ++i) {
    A.mtxIndL[i] = A.mtxIndL[0] + i * numberOfNonzerosPerRow;
    A.matrixValues[i] = A.matrixValues[0] + i * numberOfNonzerosPerRow;
    A.mtxIndG[i] = A.mtxIndG[0] + i * numberOfNonzerosPerRow;
  }

  ifile.read(reinterpret_cast<char*>(&A.mtxIndG[0][0]), A.localNumberOfRows*numberOfNonzerosPerRow*sizeof(A.mtxIndG[0][0]));
  ifile.read(reinterpret_cast<char*>(&A.mtxIndL[0][0]), A.localNumberOfRows*numberOfNonzerosPerRow*sizeof(A.mtxIndL[0][0]));
  ifile.read(reinterpret_cast<char*>(&A.matrixValues[0][0]), A.localNumberOfRows*numberOfNonzerosPerRow*sizeof(A.matrixValues[0][0]));
#endif
  // matrixDiagonal points to matrixValues, read offsets
  // A.matrixDiagonal[i] = &A.matrixValues[i][0] - offset;
  for (local_int_t i = 0; i < A.localNumberOfRows; i++) {
    ptrdiff_t offset = 0;
    ifile.read(reinterpret_cast<char*>(&offset), sizeof(offset));
    A.matrixDiagonal[i] = &A.matrixValues[i][0] - offset;
  }
  // for (int i = 0; i < A.geom->nx*A.geom->ny*A.geom->nz; ++i) {
  //   global_int_t key = 0;
  //   local_int_t value = 0;
  //   ifile.read(reinterpret_cast<char*>(&key), sizeof(key));
  //   ifile.read(reinterpret_cast<char*>(&value), sizeof(value));
  //   A.globalToLocalMap[key] = value;
  // }
  A.localToGlobalMap.resize(A.localNumberOfRows);
  ifile.read(reinterpret_cast<char*>(&A.localToGlobalMap[0]), A.localNumberOfRows*sizeof(A.localToGlobalMap[0]));
  ifile.read(reinterpret_cast<char*>(&A.isDotProductOptimized), sizeof(A.isDotProductOptimized));
  ifile.read(reinterpret_cast<char*>(&A.isSpmvOptimized), sizeof(A.isSpmvOptimized));
  ifile.read(reinterpret_cast<char*>(&A.isMgOptimized), sizeof(A.isMgOptimized));
  ifile.read(reinterpret_cast<char*>(&A.isWaxpbyOptimized), sizeof(A.isWaxpbyOptimized));
  ifile.read(reinterpret_cast<char*>(&A.numberOfExternalValues), sizeof(A.numberOfExternalValues));
  ifile.read(reinterpret_cast<char*>(&A.numberOfSendNeighbors), sizeof(A.numberOfSendNeighbors));
  ifile.read(reinterpret_cast<char*>(&A.totalToBeSent), sizeof(A.totalToBeSent));
  A.elementsToSend = new local_int_t[A.totalToBeSent];
  ifile.read(reinterpret_cast<char*>(&A.elementsToSend[0]), A.totalToBeSent*sizeof(A.elementsToSend[0]));
  A.neighbors = new int[A.numberOfSendNeighbors];
  ifile.read(reinterpret_cast<char*>(&A.neighbors[0]), A.numberOfSendNeighbors*sizeof(A.neighbors[0]));
  A.receiveLength = new local_int_t[A.numberOfSendNeighbors];
  ifile.read(reinterpret_cast<char*>(&A.receiveLength[0]), A.numberOfSendNeighbors*sizeof(A.receiveLength[0]));
  A.sendLength = new local_int_t[A.numberOfSendNeighbors];
  ifile.read(reinterpret_cast<char*>(&A.sendLength[0]), A.numberOfSendNeighbors*sizeof(A.sendLength[0]));
  A.sendBuffer = new double[A.totalToBeSent];
  ifile.close();
}

static void VectorSerializeImpl(Vector &v, std::ofstream &ofile) {
  ofile.write(reinterpret_cast<char*>(&v.localLength), sizeof(v.localLength));
  ofile.write(reinterpret_cast<char*>(&v.values[0]), v.localLength*sizeof(v.values[0]));
}

static void VectorDeserializeImpl(Vector &v, std::ifstream &ifile) {
  ifile.read(reinterpret_cast<char*>(&v.localLength), sizeof(v.localLength));
  InitializeVector(v, v.localLength);
  ifile.read(reinterpret_cast<char*>(&v.values[0]), v.localLength*sizeof(v.values[0]));
}

void MGDataSerialize(MGData &mgData, const std::string &fname) {
  std::ofstream ofile;
  ofile.open(fname.c_str(), std::ios::out | std::ios::binary | std::ios::trunc);
  die_if(ofile.fail(), "MGData: open() error");
  ofile.write(reinterpret_cast<char*>(&mgData.numberOfPresmootherSteps), sizeof(mgData.numberOfPresmootherSteps));
  ofile.write(reinterpret_cast<char*>(&mgData.numberOfPostsmootherSteps), sizeof(mgData.numberOfPostsmootherSteps));

  ofile.write(reinterpret_cast<char*>(&mgData.f2cOperatorLength), sizeof(mgData.f2cOperatorLength));
  ofile.write(reinterpret_cast<char*>(&mgData.f2cOperator[0]), mgData.f2cOperatorLength*sizeof(mgData.f2cOperator[0]));

  VectorSerializeImpl(*mgData.rc, ofile);
  VectorSerializeImpl(*mgData.xc, ofile);
  VectorSerializeImpl(*mgData.Axf, ofile);

  ofile.close();
}

void MGDataDeserialize(MGData &mgData, const std::string &fname) {
  std::ifstream ifile;
  ifile.open(fname.c_str(), std::ios::in | std::ios::binary);
  die_if(ifile.fail(), "MGData: open() error");
  ifile.read(reinterpret_cast<char*>(&mgData.numberOfPresmootherSteps), sizeof(mgData.numberOfPresmootherSteps));
  ifile.read(reinterpret_cast<char*>(&mgData.numberOfPostsmootherSteps), sizeof(mgData.numberOfPostsmootherSteps));
  ifile.read(reinterpret_cast<char*>(&mgData.f2cOperatorLength), sizeof(mgData.f2cOperatorLength));
  mgData.f2cOperator = new local_int_t[mgData.f2cOperatorLength];
  ifile.read(reinterpret_cast<char*>(&mgData.f2cOperator[0]), mgData.f2cOperatorLength*sizeof(mgData.f2cOperator[0]));

  mgData.rc = new Vector;
  VectorDeserializeImpl(*mgData.rc, ifile);

  mgData.xc = new Vector;
  VectorDeserializeImpl(*mgData.xc, ifile);

  mgData.Axf = new Vector;
  VectorDeserializeImpl(*mgData.Axf, ifile);

  ifile.close();
}

void VectorSerialize(Vector &v, const std::string &fname) {
  std::ofstream ofile;
  ofile.open(fname.c_str(), std::ios::out | std::ios::binary | std::ios::trunc);
  die_if(ofile.fail(), "Vector: open() error");

  VectorSerializeImpl(v, ofile);

  ofile.close();
}

void VectorDeserialize(Vector &v, const std::string &fname) {
  std::ifstream ifile;
  ifile.open(fname.c_str(), std::ios::in | std::ios::binary);
  die_if(ifile.fail(), "Vector: open() error");

  VectorDeserializeImpl(v, ifile);

  ifile.close();
}

void ToleranceSerialize(double &d, const std::string &fname) {
  std::ofstream ofile;
  ofile.open(fname.c_str(), std::ios::out | std::ios::binary | std::ios::trunc);
  die_if(ofile.fail(), "Tolerance: open() error");
  ofile.write(reinterpret_cast<char*>(&d), sizeof(d));
  ofile.close();
}

void ToleranceDeserialize(double &d, const std::string &fname) {
  std::ifstream ifile;
  ifile.open(fname.c_str(), std::ios::in | std::ios::binary);
  die_if(ifile.fail(), "Tolerance: open() error");

  ifile.read(reinterpret_cast<char*>(&d), sizeof(d));

  ifile.close();
}
