/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://github.com/gstlearn                                       */
/* License: BSD 3 clauses                                                     */
/*                                                                            */
/******************************************************************************/
#pragma once

#include "Basic/VectorNumT.hpp"

#include <Eigen/Dense>

#include "gstlearn_export.hpp"

class GSTLEARN_EXPORT MatrixEigen
{
public:

	MatrixEigen(int i = 1,int j = 1,double val = 0.);
	MatrixEigen(const VectorVectorDouble&,bool transpose = false);
	MatrixEigen(const Eigen::MatrixXd&);
	MatrixEigen sumRows() const;
	MatrixEigen sumCols() const;
	void sumColsInPlace(VectorDouble& res) const;
	void sumRowsInPlace(VectorDouble& res) const;
	const Eigen::MatrixXd* getMatrix() const {return &_matrix;}
	Eigen::MatrixXd* getMatrixM() {return &_matrix;}
	double* data(){return _matrix.data();}
	void setCol(int j,const VectorDouble& vect);
	void setRow(int j,const VectorDouble& vect);
	VectorDouble convert2VectorDouble() const;
	VectorDouble prodMatVec(VectorDouble& in)const;
	VectorDouble prodTMatVec(VectorDouble& in)const;
	void prodMatVecInPlace(VectorDouble& in,VectorDouble& out)const;
	void prodTMatVecInPlace(VectorDouble& in,VectorDouble& out)const;
	void prodTMatVecInPlace(const MatrixEigen &in,MatrixEigen& out)const;

	static MatrixEigen prod(const MatrixEigen& mat1,const MatrixEigen& mat2);
	static MatrixEigen prodT1(const MatrixEigen& mat1,const MatrixEigen& mat2);
	static MatrixEigen prodT2(const MatrixEigen& mat1,const MatrixEigen& mat2);
	static MatrixEigen productPointwise(const MatrixEigen& mat1,const MatrixEigen& mat2);
	double get(int i, int j)const {return _matrix(i,j);}
	void set(int i, int j, double val){ _matrix(i,j)=val;}
	MatrixEigen solve(MatrixEigen& rhs) const;
	void solveByChol(const VectorDouble& rhs,VectorDouble& res) const;
	void solveByChol(const MatrixEigen& rhs,MatrixEigen& res) const;
	void sumElem(int,int,double);
	int rows() const {return _matrix.rows();}
	int cols() const {return _matrix.cols();}
    virtual ~MatrixEigen();
    void computeInverse() const;


    void solveInPlace(const MatrixEigen& in,MatrixEigen& res) const;

private:
    void _prepareFactor() const;
    void _solveByChol(const Eigen::VectorXd& rhs,Eigen::VectorXd& res) const;
    void _solveByChol(const Eigen::MatrixXd& rhs,Eigen::MatrixXd& res) const;

    void _prepareInverse()const;
    void _reset();


private:
    Eigen::MatrixXd _matrix;
    mutable Eigen::LLT<Eigen::MatrixXd> _factor;
    mutable Eigen::MatrixXd _inverse;
    mutable bool _factorComputed;
    mutable bool _inverseComputed;
};