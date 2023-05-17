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
MatrixEigen::MatrixEigen(int n,int p,double val) {

	 omp_set_num_threads(1);
	  _matrix= Eigen::MatrixXd::Constant(n,p,val);

}

MatrixEigen::MatrixEigen(const VectorVectorDouble& mat,bool transpose)
{


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
	_matrix = mat;
}
void MatrixEigen::setCol(int j,const VectorDouble& vect)
{
	for (int i = 0; i < (int)vect.size();i++)
	{
		_matrix(i,j) = vect[i];
	}
}

void MatrixEigen::sumElem(int i,int j,double val)
{
	_matrix(i,j) += val;
}

void MatrixEigen::prodMatVecInPlace(VectorDouble& in,VectorDouble& out)
{

	Eigen::Map<Eigen::VectorXd> inv(in.data(), in.size());
	Eigen::Map<Eigen::VectorXd> outv(out.data(), out.size());

	outv.noalias() += _matrix * inv;
}
void MatrixEigen::prodTMatVecInPlace(VectorDouble& in,VectorDouble& out)
{
	Eigen::Map<Eigen::VectorXd> inv(in.data(), in.size());
	Eigen::Map<Eigen::VectorXd> outv(out.data(), out.size());

	outv.noalias() += _matrix.transpose() * inv;
}

MatrixEigen MatrixEigen::solve(MatrixEigen& rhs) const
{

	MatrixEigen result = MatrixEigen(_matrix.rows(),_matrix.cols());
	result._matrix = _matrix.inverse() * rhs._matrix;
	return result;
}

void MatrixEigen::solve(VectorDouble& rhs,VectorDouble& res) const
{
	Eigen::Map<Eigen::VectorXd> rhsv(rhs.data(), rhs.size());
	Eigen::Map<Eigen::VectorXd> resv(res.data(), res.size());
	resv.noalias() += _matrix.llt().solve(rhsv);

}

void MatrixEigen::setRow(int j,const VectorDouble& vect)
{
	for (int i = 0; i < (int)vect.size();i++)
	{
		_matrix(j,i) = vect[i];
	}
}


VectorDouble MatrixEigen::prodMatVec(VectorDouble& in)
{
	VectorDouble result = VectorDouble(_matrix.rows());
	prodMatVecInPlace(in, result);
	return result;
}

VectorDouble MatrixEigen::prodTMatVec(VectorDouble& in)
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

MatrixEigen::~MatrixEigen() {
	// TODO Auto-generated destructor stub
}

