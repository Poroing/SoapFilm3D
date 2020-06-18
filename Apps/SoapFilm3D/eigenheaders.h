//
//  eigenheaders.h
//  MultiTracker
//
//  Created by Fang Da on 10/24/14.
//
//

#ifndef MultiTracker_eigenheaders_h
#define MultiTracker_eigenheaders_h

#include "surftrack.h"
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <boost/range/irange.hpp>

typedef Eigen::Matrix<double, 4, 4> Mat4d;
typedef Eigen::Matrix<double, 3, 3> Mat3d;
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> MatXd;
typedef Eigen::Matrix<double, 4, 1> Vec4d;
typedef Eigen::Matrix<double, 3, 1> Vec3d;
typedef Eigen::Matrix<double, 2, 1> Vec2d;
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> VecXd;

typedef Eigen::Matrix<int, 3, 1> Vec3i;
typedef Eigen::Matrix<int, 2, 1> Vec2i;
typedef Eigen::Matrix<int, Eigen::Dynamic, 3> MatXi;

typedef Eigen::SparseMatrix<double> SparseMatd;

Vec3i
vc(const LosTopos::Vec3st& t);
Vec2i
vc(const LosTopos::Vec2i& l);
Vec3d
vc(const LosTopos::Vec3d& v);
LosTopos::Vec3d
vc(const Vec3d& v);

template<typename Derived>
class MatrixComp
{
  public:
    bool operator()(const Eigen::PlainObjectBase<Derived>& lhs,
                    const Eigen::PlainObjectBase<Derived>& rhs) const
    {
        assert(lhs.cols() == rhs.cols());
        assert(lhs.rows() == rhs.rows());

        for (size_t col : boost::irange((typename Derived::Index)0, lhs.cols()))
        {
            for (size_t row : boost::irange((typename Derived::Index)0, rhs.rows()))
            {
                if (lhs(row, col) < rhs(row, col))
                {
                    return true;
                }
                if (lhs(row, col) > rhs(row, col))
                {
                    return false;
                }
            }
        }

        return false;
    }
};

using Vec2iComp = MatrixComp<Vec2i>;

#endif
