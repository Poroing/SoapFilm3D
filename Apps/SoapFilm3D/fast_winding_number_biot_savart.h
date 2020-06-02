#ifndef IGL_FAST_WINDING_NUMBER
#define IGL_FAST_WINDING_NUMBER
#include "igl/igl_inline.h"
#include <Eigen/Core>
#include <vector>
#include <Eigen/Geometry>
#include "igl/octree.h"
#include "igl/knn.h"
#include "igl/parallel_for.h"
#include "igl/PI.h"
#include "igl/cross.h"
#include <vector>
#include <cassert>

namespace igl
{
  // Generate the precomputation for the fast winding number for point data
  // [Barill et. al 2018].
  //
  // Given a set of 3D points P, with normals N, areas A, along with octree
  // data, and an expansion order, we define a taylor series expansion at each
  // octree cell.
  //
  // The octree data is designed to come from igl::octree, and the areas (if not
  // obtained at scan time), may be calculated using
  // igl::copyleft::cgal::point_areas.
  //
  // Inputs:
  //   P  #P by 3 list of point locations
  //   G  #P by 3 list of sheet strengths
  //   point_indices  a vector of vectors, where the ith entry is a vector of
  //                  the indices into P that are the ith octree cell's points
  //   CH             #OctreeCells by 8, where the ith row is the indices of
  //                  the ith octree cell's children
  //   expansion_order    the order of the taylor expansion. We support 0,1,2.
  // Outputs:
  //   CM  #OctreeCells by 3 list of each cell's center of mass
  //   R   #OctreeCells by 1 list of each cell's maximum distance of any point
  //       to the center of mass
  //   EC  #OctreeCells by #TaylorCoefficients list of expansion coefficients.
  //       (Note that #TaylorCoefficients = ∑_{i=1}^{expansion_order} 3^i)
  template <
    typename DerivedP, 
    typename DerivedG,
    typename Index, 
    typename DerivedCH, 
    typename DerivedCM, 
    typename DerivedR,
    typename DerivedEC>
  IGL_INLINE void fast_winding_number_biot_savart(
    const Eigen::MatrixBase<DerivedP>& P,
    const Eigen::MatrixBase<DerivedG>& G,
    const std::vector<std::vector<Index> > & point_indices,
    const Eigen::MatrixBase<DerivedCH>& CH,
    int expansion_order,
    Eigen::PlainObjectBase<DerivedCM>& CM,
    Eigen::PlainObjectBase<DerivedR>& R,
    Eigen::PlainObjectBase<DerivedEC>& EC);
  // Evaluate the fast winding number for point data, having already done the
  // the precomputation
  //
  // Inputs:
  //   P  #P by 3 list of point locations
  //   point_indices  a vector of vectors, where the ith entry is a vector of
  //                  the indices into P that are the ith octree cell's points
  //   CH  #OctreeCells by 8, where the ith row is the indices of
  //       the ith octree cell's children
  //   CM  #OctreeCells by 3 list of each cell's center of mass
  //   R   #OctreeCells by 1 list of each cell's maximum distance of any point
  //       to the center of mass
  //   EC  #OctreeCells by #TaylorCoefficients list of expansion coefficients.
  //        (Note that #TaylorCoefficients = ∑_{i=1}^{expansion_order} 3^i)
  //   Q  #Q by 3 list of query points for the winding number
  //   beta  This is a Barnes-Hut style accuracy term that separates near feild
  //         from far field. The higher the beta, the more accurate and slower
  //         the evaluation. We reccommend using a beta value of 2. Note that
  //         for a beta value ≤ 0, we use the direct evaluation, rather than
  //         the fast approximation
  // Outputs:
  //   WN  #Q by 1 list of windinng number values at each query point
  //
  template <
    typename DerivedP, 
    typename DerivedG,
    typename Index, 
    typename DerivedCH, 
    typename DerivedCM, 
    typename DerivedR,
    typename DerivedEC, 
    typename DerivedQ, 
    typename BetaType,
    typename DeltaType,
    typename DerivedWN>
  IGL_INLINE void fast_winding_number_biot_savart(
    const Eigen::MatrixBase<DerivedP>& P,
    const Eigen::MatrixBase<DerivedG>& G,
    const std::vector<std::vector<Index> > & point_indices,
    const Eigen::MatrixBase<DerivedCH>& CH,
    const Eigen::MatrixBase<DerivedCM>& CM,
    const Eigen::MatrixBase<DerivedR>& R,
    const Eigen::MatrixBase<DerivedEC>& EC,
    const Eigen::MatrixBase<DerivedQ>& Q,
    const BetaType beta,
    DeltaType delta,
    Eigen::PlainObjectBase<DerivedWN>& WN);
  template <
    typename DerivedP, 
    typename DerivedG,
    typename DerivedQ, 
    typename BetaType, 
    typename DeltaType,
    typename DerivedWN>
  IGL_INLINE void fast_winding_number_biot_savart(
    const Eigen::MatrixBase<DerivedP>& P,
    const Eigen::MatrixBase<DerivedG>& G,
    const Eigen::MatrixBase<DerivedQ>& Q,
    const int expansion_order,
    const BetaType beta,
    DeltaType delta,
    Eigen::PlainObjectBase<DerivedWN>& WN);
  // Evaluate the fast winding number for point data, with default expansion
  // order and beta (both are set to 2).
  //
  // This function performes the precomputation and evaluation all in one.
  // If you need to acess the precomuptation for repeated evaluations, use the
  // two functions designed for exposed precomputation (described above).
  //
  // Inputs:
  //   P  #P by 3 list of point locations
  //   Q  #Q by 3 list of query points for the winding number
  // Outputs:
  //   WN  #Q by 1 list of windinng number values at each query point
  //
  template <
    typename DerivedP, 
    typename DerivedG,
    typename DerivedQ, 
    typename DerivedWN,
    typename DeltaType>
  IGL_INLINE void fast_winding_number_biot_savart(
    const Eigen::MatrixBase<DerivedP>& P,
    const Eigen::MatrixBase<DerivedG>& G,
    const Eigen::MatrixBase<DerivedQ>& Q,
    DeltaType delta,
    Eigen::PlainObjectBase<DerivedWN>& WN);
}

template <
  typename DerivedP, 
  typename DerivedG,
  typename Index, 
  typename DerivedCH, 
  typename DerivedCM, 
  typename DerivedR,
  typename DerivedEC>
IGL_INLINE void igl::fast_winding_number_biot_savart(
  const Eigen::MatrixBase<DerivedP>& P,
  const Eigen::MatrixBase<DerivedG>& G,
  const std::vector<std::vector<Index> > & point_indices,
  const Eigen::MatrixBase<DerivedCH>& CH,
  int expansion_order,
  Eigen::PlainObjectBase<DerivedCM>& CM,
  Eigen::PlainObjectBase<DerivedR>& R,
  Eigen::PlainObjectBase<DerivedEC>& EC)
{
  typedef typename DerivedP::Scalar real_p;
  typedef typename DerivedCM::Scalar real_cm;
  typedef typename DerivedR::Scalar real_r;
  typedef typename DerivedEC::Scalar real_ec;

  typedef Eigen::Matrix<real_p,1,3> RowVec3p;

  int m = CH.size();
  int num_terms;

  assert(expansion_order < 3 && expansion_order >= 0 && "m must be less than n");
  if(expansion_order == 0){
      num_terms = 3;
  } else if(expansion_order ==1){
      num_terms = 3 + 9;
  } else if(expansion_order == 2){
      num_terms = 3 + 9 + 27;
  }

  R.resize(m);
  CM.resize(m,3);
  EC.resize(m,num_terms);
  EC.setZero(m,num_terms);
  std::function< void(const int) > helper;
  helper = [&helper,
            &P,&G,&expansion_order,&point_indices,&CH,&EC,&R,&CM]
  (const int index)-> void
  {
      Eigen::Matrix<real_cm,1,3> masscenter;
      masscenter << 0,0,0;
      Eigen::Matrix<real_ec,1,3> zeroth_expansion;
      zeroth_expansion << 0,0,0;
      real_p areatotal = 0.0;

      for(int j = 0; j < point_indices[index].size(); j++){
          int curr_point_index = point_indices[index][j];
        
          masscenter += P.row(curr_point_index);
          zeroth_expansion += G.row(curr_point_index);
          areatotal += 1.;
      }
    
      masscenter = masscenter / areatotal;
      CM.row(index) = masscenter;
      EC.block(index,0,1,3) = zeroth_expansion;
    
      real_r max_norm = 0;
      real_r curr_norm;
    
      for(int i = 0; i < point_indices[index].size(); i++){
          //Get max distance from center of mass:
          int curr_point_index = point_indices[index][i];
          Eigen::Matrix<real_r,1,3> point =
              P.row(curr_point_index)-masscenter;
          curr_norm = point.norm();
          if(curr_norm > max_norm){
              max_norm = curr_norm;
          }
        
          //Calculate higher order terms if necessary
          Eigen::Matrix<real_ec,3,3> TempCoeffs;
          if(EC.cols() >= (3+9)){
              TempCoeffs = point.transpose()*
                              G.row(curr_point_index);
              EC.block(index,3,1,9) +=
              Eigen::Map<Eigen::Matrix<real_ec,1,9> >(TempCoeffs.data(),
                                                      TempCoeffs.size());
          }
        
          if(EC.cols() == (3+9+27)){
              for(int k = 0; k < 3; k++){
                  TempCoeffs = 0.5 * point(k) * (point.transpose()*G.row(curr_point_index));
                  EC.block(index,12+9*k,1,9) += Eigen::Map<
                    Eigen::Matrix<real_ec,1,9> >(TempCoeffs.data(),
                                                 TempCoeffs.size());
              }
          }
      }
    
      R(index) = max_norm;
      if(CH(index,0) != -1)
      {
          for(int i = 0; i < 8; i++){
              int child = CH(index,i);
              helper(child);
          }
      }
  };
  helper(0);
}

template <
  typename DerivedP, 
  typename DerivedG,
  typename Index, 
  typename DerivedCH, 
  typename DerivedCM, 
  typename DerivedR,
  typename DerivedEC, 
  typename DerivedQ, 
  typename BetaType,
  typename DeltaType,
  typename DerivedWN>
IGL_INLINE void igl::fast_winding_number_biot_savart(
  const Eigen::MatrixBase<DerivedP>& P,
  const Eigen::MatrixBase<DerivedG>& G,
  const std::vector<std::vector<Index> > & point_indices,
  const Eigen::MatrixBase<DerivedCH>& CH,
  const Eigen::MatrixBase<DerivedCM>& CM,
  const Eigen::MatrixBase<DerivedR>& R,
  const Eigen::MatrixBase<DerivedEC>& EC,
  const Eigen::MatrixBase<DerivedQ>& Q,
  const BetaType beta,
  DeltaType delta,
  Eigen::PlainObjectBase<DerivedWN>& WN)
{

  typedef typename DerivedP::Scalar real_p;
  typedef typename DerivedG::Scalar real_g;
  typedef typename DerivedCM::Scalar real_cm;
  typedef typename DerivedR::Scalar real_r;
  typedef typename DerivedEC::Scalar real_ec;
  typedef typename DerivedQ::Scalar real_q;
  typedef typename DerivedWN::Scalar real_wn;
  const real_wn PI_4 = 4.0*igl::PI;

  typedef Eigen::Matrix<
    typename DerivedEC::Scalar,
    1,
    DerivedEC::ColsAtCompileTime> ECRow;

  typedef Eigen::Matrix<real_q,1,3> RowVec;
  typedef Eigen::Matrix<real_wn, 1, 3> WnRowVec;
  typedef Eigen::Matrix<real_ec,3,3> EC_3by3;

  auto direct_eval = [&PI_4,delta](
    const RowVec & loc,
    const Eigen::Matrix<real_ec,1,3> & anorm)->WnRowVec
  {
    const typename RowVec::Scalar loc_norm = std::sqrt(loc.norm() * loc.norm() + delta * delta);
    WnRowVec result;
    igl::cross(anorm, loc, result);
    return - result / (PI_4*(loc_norm*loc_norm*loc_norm));
  };

  auto expansion_eval = 
    [&direct_eval,&EC,&PI_4,delta](
      const RowVec & loc,
      const int & child_index)->WnRowVec
  {
    WnRowVec wn;
    wn = direct_eval(loc,EC.row(child_index).template head<3>());
    real_wn r = std::sqrt(loc.norm() * loc.norm() + delta * delta);
    real_wn PI_4_r3;
    real_wn PI_4_r5;
    real_wn PI_4_r7;
    if(EC.row(child_index).size()>3)
    {
      PI_4_r3 = PI_4*r*r*r;
      PI_4_r5 = PI_4_r3*r*r;
      const real_ec d = 1.0/(PI_4_r3);
      Eigen::Matrix<real_ec,3,3> SecondDerivative = 
        loc.transpose()*loc*(-3.0/(PI_4_r5));
      SecondDerivative(0,0) += d;
      SecondDerivative(1,1) += d;
      SecondDerivative(2,2) += d;

      Eigen::Matrix<real_wn, 3, 3> TempCoeffs;
      Eigen::Matrix<real_ec, 3, 3> ExpansionCoefficients = Eigen::Map<const Eigen::Matrix<real_ec, 3, 3>>(EC.row(child_index).template segment<9>(3).data());
      igl::cross(
              ExpansionCoefficients,
              SecondDerivative,
              TempCoeffs);
      for (size_t i = 0; i < 3; ++i)
      {
          wn += TempCoeffs.row(i);
      }
    }
    if(EC.row(child_index).size()>3+9)
    {
      PI_4_r7 = PI_4_r5*r*r;
      const Eigen::Matrix<real_ec,3,3> locTloc = loc.transpose()*(loc/(PI_4_r7));
      for(int i = 0; i < 3; i++)
      {
        Eigen::Matrix<real_ec,3,3> RowCol_Diagonal = 
          Eigen::Matrix<real_ec,3,3>::Zero(3,3);
        for(int u = 0;u<3;u++)
        {
          for(int v = 0;v<3;v++)
          {
            if(u==v) RowCol_Diagonal(u,v) += loc(i);
            if(u==i) RowCol_Diagonal(u,v) += loc(v);
            if(v==i) RowCol_Diagonal(u,v) += loc(u);
          }
        }
        Eigen::Matrix<real_ec,3,3> ThirdDerivative = 
          15.0*loc(i)*locTloc + (-3.0/(PI_4_r5))*(RowCol_Diagonal);

        Eigen::Matrix<real_wn, 3, 3> TempCoeffs;
        Eigen::Matrix<real_ec, 3, 3> ExpansionCoefficients = Eigen::Map<const Eigen::Matrix<real_ec,3,3>>(EC.row(child_index).template segment<9>(12 + i*9).data());
        igl::cross(
                ExpansionCoefficients,
                ThirdDerivative,
                TempCoeffs);

        for (size_t j = 0; j < 3; ++j)
        {
            wn -= TempCoeffs.row(j);
        }
      }
    }
    return wn;
  };

  int m = Q.rows();
  WN.resize(m,3);

  std::function< WnRowVec(const RowVec & , const std::vector<int> &) > helper;
  helper = [&helper,
            &P, &G,
            &point_indices,&CH,
            &CM,&R,&EC,&beta,
            &direct_eval,&expansion_eval]
  (const RowVec & query, const std::vector<int> & near_indices)-> WnRowVec
  {
    WnRowVec wn(WnRowVec::Zero());
    std::vector<int> new_near_indices;
    new_near_indices.reserve(8);
    for(int i = 0; i < near_indices.size(); i++)
    {
      int index = near_indices[i];
      //Leaf Case, Brute force
      if(CH(index,0) == -1)
      {
        for(int j = 0; j < point_indices[index].size(); j++)
        {
          int curr_row = point_indices[index][j];
          wn += direct_eval(P.row(curr_row)-query,
                            G.row(curr_row));
        }
      }
      //Non-Leaf Case
      else 
      {
        for(int child = 0; child < 8; child++)
        {
          int child_index = CH(index,child);
          if(point_indices[child_index].size() > 0)
          {
            const RowVec CMciq = (CM.row(child_index)-query);
            if(CMciq.norm() > beta*R(child_index))
            {
              if(CH(child_index,0) == -1)
              {
                for(int j=0;j<point_indices[child_index].size();j++)
                {
                  int curr_row = point_indices[child_index][j];
                  wn += direct_eval(P.row(curr_row)-query,
                                    G.row(curr_row));
                }
              }else{
                wn += expansion_eval(CMciq,child_index);
              }
            }else 
            {
              new_near_indices.emplace_back(child_index);
            }
          }
        }
      }
    }
    if(new_near_indices.size() > 0)
    {
      wn += helper(query,new_near_indices);
    }
    return wn;
  };

  if(beta > 0)
  {
    const std::vector<int> near_indices_start = {0};
    igl::parallel_for(m,[&](int iter){
      WN.row(iter) = helper(Q.row(iter).eval(), near_indices_start);
    },1000);
  } else 
  {
    igl::parallel_for(m,[&](int iter){
      WnRowVec wn(WnRowVec::Zero());
      for(int j = 0; j < P.rows(); j++)
      {
        wn += direct_eval(P.row(j) - Q.row(iter), G.row(j));
      }
      WN.row(iter) = wn;
    },1000);
  }
}

template <
  typename DerivedP, 
  typename DerivedG, 
  typename DerivedQ, 
  typename BetaType, 
  typename DeltaType,
  typename DerivedWN>
IGL_INLINE void igl::fast_winding_number_biot_savart(
  const Eigen::MatrixBase<DerivedP>& P,
  const Eigen::MatrixBase<DerivedG>& G,
  const Eigen::MatrixBase<DerivedQ>& Q,
  const int expansion_order,
  const BetaType beta,
  DeltaType delta,
  Eigen::PlainObjectBase<DerivedWN>& WN)
{
  typedef typename DerivedWN::Scalar real;
  
  std::vector<std::vector<int> > point_indices;
  Eigen::Matrix<int,Eigen::Dynamic,8> CH;
  Eigen::Matrix<real,Eigen::Dynamic,3> CN;
  Eigen::Matrix<real,Eigen::Dynamic,1> W;

  octree(P,point_indices,CH,CN,W);

  Eigen::Matrix<real,Eigen::Dynamic,Eigen::Dynamic> EC;
  Eigen::Matrix<real,Eigen::Dynamic,3> CM;
  Eigen::Matrix<real,Eigen::Dynamic,1> R;

  fast_winding_number_biot_savart(P,G,point_indices,CH,expansion_order,CM,R,EC);
  fast_winding_number_biot_savart(P,G,point_indices,CH,CM,R,EC,Q,beta,delta,WN);
}

template <
  typename DerivedP, 
  typename DerivedG, 
  typename DerivedQ, 
  typename DerivedWN,
  typename DeltaType>
IGL_INLINE void igl::fast_winding_number_biot_savart(
  const Eigen::MatrixBase<DerivedP>& P,
  const Eigen::MatrixBase<DerivedG>& G,
  const Eigen::MatrixBase<DerivedQ>& Q,
  DeltaType delta,
  Eigen::PlainObjectBase<DerivedWN>& WN)
{
  fast_winding_number_biot_savart(P,G,Q,2,2.0,delta,WN);
}

#endif

