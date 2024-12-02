#ifndef SERIALIZATION_HPP
#define SERIALIZATION_HPP

#include <string>

void ParamsSerialize(HPCG_Params &params, const std::string &fname = "HPCG_Params.bin");
void ParamsDeserialize(HPCG_Params &params, const std::string &fname = "HPCG_Params.bin");

void GeometrySerialize(Geometry &geom, const std::string &fname = "Geometry.bin");
void GeometryDeserialize(Geometry &geom, const std::string &fname = "Geometry.bin");

void SparseMatrixSerialize(SparseMatrix &A, const std::string &fname = "SparseMatrix.bin");
void SparseMatrixDeserialize(SparseMatrix &A, const std::string &fname = "SparseMatrix.bin");

void MGDataSerialize(MGData &A, const std::string &fname = "MGData.bin");
void MGDataDeserialize(MGData &A, const std::string &fname = "MGData.bin");

void VectorSerialize(Vector &v, const std::string &fname);
void VectorDeserialize(Vector &v, const std::string &fname);

void ToleranceSerialize(double &d, const std::string &fname = "refTolerance.bin");
void ToleranceDeserialize(double &d, const std::string &fname = "refTolerance.bin");

#endif // SERIALIZATION_HPP

