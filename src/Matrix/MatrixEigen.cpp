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
#include "Matrix/MatrixEigen.hpp"
#include <Eigen/Cholesky>

#include <iostream>

MatrixEigen::MatrixEigen(int n,int p,double val)
{

	_reset();
	  _matrix= Eigen::MatrixXd::Constant(n,p,val);

}

void MatrixEigen::_reset()
{
	  _factorComputed = false;
	  _inverseComputed = false;
}

MatrixEigen::MatrixEigen(const VectorVectorDouble& mat,bool transpose)
{
	_reset();

	if (!transpose)
	{
		_matrix = Eigen::MatrixXd((int)mat.size(),(int)mat[0].size());
		for(int i = 0; i < (int)mat.size();i++)
		{
			for (int j = 0; j < (int) mat[0].size(); j++)
			{
				_matrix(i,j) = mat[i][j];
			}
		}
	}

	if (transpose)
	{
		_matrix = Eigen::MatrixXd((int)mat[0].size(),(int)mat.size());
		for(int i = 0; i < (int)mat.size();i++)
		{
			for (int j = 0; j < (int) mat[0].size(); j++)
			{
				_matrix(j,i) = mat[i][j];
			}
		}
	}
}

MatrixEigen::MatrixEigen(const Eigen::MatrixXd& mat)
{
	_reset();
	_matrix = mat;
}


MatrixEigen MatrixEigen::sumRows() const
{
	return MatrixEigen(_matrix.rowwise().sum());

}
MatrixEigen MatrixEigen::sumCols() const
{
	return MatrixEigen(_matrix.colwise().sum());

}

void MatrixEigen::sumColsInPlace(VectorDouble& res) const
{

	Eigen::Map<Eigen::VectorXd> result(res.data(), res.size());
	result.noalias() += _matrix.colwise().sum();

}

void MatrixEigen::sumRowsInPlace(VectorDouble& res) const
{

	Eigen::Map<Eigen::VectorXd> result(res.data(), res.size());
	result.noalias() += _matrix.rowwise().sum();

}


void MatrixEigen::setCol(int j,const VectorDouble& vect)
{
	_reset();
	for (int i = 0; i < (int)vect.size();i++)
	{
		_matrix(i,j) = vect[i];
	}
}

void MatrixEigen::sumElem(int i,int j,double val)
{
	_reset();
	_matrix(i,j) += val;
}

void MatrixEigen::prodMatVecInPlace(VectorDouble& in,VectorDouble& out)const
{

	Eigen::Map<Eigen::VectorXd> inv(in.data(), in.size());
	Eigen::Map<Eigen::VectorXd> outv(out.data(), out.size());

	outv.noalias() += _matrix * inv;
}
void MatrixEigen::prodTMatVecInPlace(VectorDouble& in,VectorDouble& out)const
{
	Eigen::Map<Eigen::VectorXd> inv(in.data(), in.size());
	Eigen::Map<Eigen::VectorXd> outv(out.data(), out.size());

	outv.noalias() += _matrix.transpose() * inv;
}

MatrixEigen MatrixEigen::solve(MatrixEigen& rhs) const
{

	MatrixEigen result = MatrixEigen(_matrix.rows(),rhs._matrix.cols());
	_prepareInverse();
	result._matrix.noalias() += _inverse * rhs._matrix;
	return result;
}

void MatrixEigen::solve(const VectorDouble& rhs,VectorDouble& res) const
{
	Eigen::Map<const Eigen::VectorXd> rhsv(rhs.data(), rhs.size());
	Eigen::Map<Eigen::VectorXd> resv(res.data(), res.size());
	_prepareFactor();
	resv.noalias() += _factor.solve(rhsv);

}

void MatrixEigen::setRow(int j,const VectorDouble& vect)
{
	_reset();
	for (int i = 0; i < (int)vect.size();i++)
	{
		_matrix(j,i) = vect[i];
	}
}


VectorDouble MatrixEigen::prodMatVec(VectorDouble& in)const
{
	VectorDouble result = VectorDouble(_matrix.rows());
	prodMatVecInPlace(in, result);
	return result;
}

VectorDouble MatrixEigen::prodTMatVec(VectorDouble& in)const
{
	VectorDouble result = VectorDouble(_matrix.rows());
	prodTMatVecInPlace(in, result);
	return result;
}
VectorDouble MatrixEigen::convert2VectorDouble()const
{
	int n = _matrix.rows() * _matrix.cols();
	std::cout<<n <<std::endl;
	VectorDouble result(_matrix.data(), _matrix.data()+n);
	return result;
}

void MatrixEigen::_prepareFactor() const
{
	if(!_factorComputed)
	{
		_factor = _matrix.llt();
		_factorComputed = true;
	}
}

void MatrixEigen::_prepareInverse() const
{
	if(!_inverseComputed)
	{
		_inverse = _matrix.inverse();
		_inverseComputed = true;
	}
}

///////////////////// static functions /////////////////////////////

MatrixEigen MatrixEigen::prod(const MatrixEigen& mat1,const MatrixEigen& mat2)
{
	MatrixEigen  res = MatrixEigen(mat1._matrix * mat2._matrix);
	return res;
}

MatrixEigen MatrixEigen::prodT1(const MatrixEigen& mat1,const MatrixEigen& mat2)
{
	MatrixEigen  res = MatrixEigen(mat1._matrix.transpose() * mat2._matrix);
	return res;
}


MatrixEigen MatrixEigen::prodT2(const MatrixEigen& mat1,const MatrixEigen& mat2)
{
	MatrixEigen  res = MatrixEigen(mat1._matrix * mat2._matrix.transpose());
	return res;
}

MatrixEigen MatrixEigen::productPointwise(const MatrixEigen& mat1,const MatrixEigen& mat2)
{
	MatrixEigen  res = MatrixEigen(mat1.rows(),mat1.cols());
	res._matrix.array() = mat1._matrix.array() * mat2._matrix.array();
	return res;
}



MatrixEigen::~MatrixEigen() {
	// TODO Auto-generated destructor stub
}

