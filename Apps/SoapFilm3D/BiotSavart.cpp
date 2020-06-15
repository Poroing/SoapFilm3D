#include "fmmtl/fmmtl/util/Clock.hpp"
#include "fast_winding_number_biot_savart.h"
#include "SimOptions.h"

#ifndef WIN32
#include "fmmtl/fmmtl/Direct.hpp"
#include "fmmtl/fmmtl/KernelMatrix.hpp"
#include "fmmtl/kernel/BiotSpherical.hpp"
#include "fmmtl/kernel/RMSpherical.hpp"
#endif

#include <boost/range/algorithm/transform.hpp>
#include <boost/range/algorithm/unique.hpp>
#include <boost/range/algorithm/sort.hpp>
#include <boost/range/algorithm_ext/erase.hpp>

#include "BiotSavart.h"

VecXd
BiotSavart_naive(const std::vector<Vec3d>& sources,
                 const std::vector<Vec3d>& targets,
                 const std::vector<Vec3d>& charges,
                 double delta)
{
    assert(charges.size() == sources.size());
    VecXd vel = VecXd::Zero(targets.size() * 3);
#pragma omp parallel
    {
#pragma omp for nowait
        for (size_t target_index = 0; target_index < targets.size(); ++target_index)
        {
            Vec3d v(0, 0, 0);
            const Vec3d& target = targets[target_index];

            for (size_t source_index : boost::irange(0lu, sources.size()))
            {

                const Vec3d& source = sources[source_index];

                Vec3d dx = target - source;
                double dxn = sqrt(dx.squaredNorm() + delta * delta);

                v += charges[source_index].cross(dx) / (dxn * dxn * dxn);
            }

            v /= (4 * M_PI);
            vel.segment<3>(target_index * 3) = v;
        }
    }

    return vel;
}

#ifndef WIN32
VecXd
BiotSavart_fmmtl(const std::vector<Vec3d>& sources,
                 const std::vector<Vec3d>& targets,
                 const std::vector<Vec3d>& charges,
                 double delta)
{
    assert(charges.size() == sources.size());
    // code adapted from FMMTL example test "error_biot.cpp"

    // Init the FMM Kernel and options
    //        typedef BiotSpherical kernel_type;
    typedef RMSpherical kernel_type;

    // Init kernel
    kernel_type K(delta, Options::intValue("fmmtl-expansion-order"));

    typedef kernel_type::point_type point_type;
    typedef kernel_type::source_type source_type;
    typedef kernel_type::target_type target_type;
    typedef kernel_type::charge_type charge_type;
    typedef kernel_type::result_type result_type;

    // Init points and charges
    std::vector<source_type> fmmtl_sources(sources.size());
    std::vector<target_type> fmmtl_targets(targets.size());
    std::vector<charge_type> fmmtl_charges(charges.size());

    auto convert_to_fmmtl = [](const Vec3d& v) { return Vec<3, double>(v[0], v[1], v[2]); };
#pragma omp parallel
    {
#pragma omp for nowait
        for (size_t target_index = 0; target_index < targets.size(); ++target_index)
        {
            fmmtl_targets[target_index] = convert_to_fmmtl(targets[target_index]);
        }

        for (size_t source_index = 0; source_index < sources.size(); ++source_index)
        {
            fmmtl_sources[source_index] = convert_to_fmmtl(sources[source_index]);
            fmmtl_charges[source_index] = convert_to_fmmtl(charges[source_index]);
        }
    }

    // Build the FMM
    fmmtl::kernel_matrix<kernel_type> A = K(fmmtl_targets, fmmtl_sources);
    FMMOptions opts;
    opts.theta = Options::boolValue("fmmtl-theta");
    A.set_options(opts);

    // Execute the FMM
    std::vector<result_type> result = A * fmmtl_charges;

    VecXd vel = VecXd::Zero(targets.size() * 3);
    for (size_t target_index : boost::irange(0lu, targets.size()))
    {
        vel[target_index * 3 + 0] = result[target_index][0];
        vel[target_index * 3 + 1] = result[target_index][1];
        vel[target_index * 3 + 2] = result[target_index][2];
    }

    vel /= (4 * M_PI);

    return vel;
}
#endif

VecXd
BiotSavart_fast_winding_number(const std::vector<Vec3d>& sources,
                 const std::vector<Vec3d>& targets,
                 const std::vector<Vec3d>& charges,
                 double delta)
{
    assert(charges.size() == sources.size());

    Eigen::MatrixX3d libigl_sources(Eigen::MatrixX3d::Zero(sources.size(), 3));
    Eigen::MatrixX3d libigl_charges(Eigen::MatrixX3d::Zero(charges.size(), 3));
    Eigen::MatrixX3d libigl_targets(Eigen::MatrixX3d::Zero(targets.size(), 3));

#pragma omp parallel
    {
#pragma omp for nowait
        for (size_t target_index = 0; target_index < targets.size(); ++target_index)
        {
            libigl_targets.row(target_index) = targets[target_index];
        }

        for (size_t source_index = 0; source_index < sources.size(); ++source_index)
        {
            libigl_sources.row(source_index) = sources[source_index];
            libigl_charges.row(source_index) = charges[source_index];
        }
    }

    // Build the FMM
    Eigen::MatrixX3d result;
    igl::fast_winding_number_biot_savart(libigl_sources,
                                         libigl_charges,
                                         libigl_targets,
                                         Options::intValue("winding-expansion-order"),
                                         Options::doubleValue("winding-beta"),
                                         delta,
                                         result);

    VecXd vel = VecXd::Zero(targets.size() * 3);
    for (size_t target_index : boost::irange(0lu, targets.size()))
    {
        vel.segment<3>(target_index * 3) = result.row(target_index);
    }

    return vel;
}

VecXd
BiotSavart(VS3D& vs, const VecXd& dx)
{
    std::cout << "NumberVertices " << vs.mesh().nv() << std::endl;
    std::cout << "BoundingBoxVolume " << vs.getBoundingBoxVolume() << std::endl;
    std::cout << "NumberFaces " << vs.mesh().nt() << std::endl;


    Clock biot_savart_duration;
    // Init points and charges
    std::vector<Vec3d> sources(vs.mesh().nt());
    std::vector<Vec3d> targets(vs.mesh().nv());
    std::vector<Vec3d> charges(vs.mesh().nt());

#pragma omp parallel
    {
#pragma omp for nowait
        for (size_t vertex_index = 0; vertex_index < vs.mesh().nv(); ++vertex_index)
        {
            targets[vertex_index] = vs.pos(vertex_index);
        }

#pragma omp for nowait
        for (size_t triangle_index = 0; triangle_index < vs.mesh().nt(); ++triangle_index)
        {
            if (vs.surfTrack()->triangle_is_all_solid(triangle_index))
                continue; // all-solid faces don't contribute vorticity.

            sources[triangle_index] = vs.getTranslatedTriangleCenter(triangle_index, dx);
            charges[triangle_index] = vs.getTranslatedTriangleSheetStrength(triangle_index, dx);
        }
    }


    // open boundary extra face contributions
    if (vs.m_obefv.size() == vs.m_obefe.size() && vs.m_obefv.size() == vs.m_obefc.size())
    {
        for (size_t open_boundary_face_index : boost::irange(0lu, vs.m_obefv.size()))
        {
            charges.push_back(vs.m_obefe[open_boundary_face_index]
                              * vs.m_obefv[open_boundary_face_index]);
            sources.push_back(vs.m_obefc[open_boundary_face_index]);
        }
    }

    VecXd result;

#ifndef WIN32
    if (Options::strValue("fast-summation") == "fmmtl")
    {
        result = BiotSavart_fmmtl(sources, targets, charges, vs.delta());
    }
    else
#endif
      if (Options::strValue("fast-summation") == "naive")
    {
        result = BiotSavart_naive(sources, targets, charges, vs.delta());
    }
    else if (Options::strValue("fast-summation") == "winding")
    {
        result = BiotSavart_fast_winding_number(sources, targets, charges, vs.delta());
    }
    else {
        throw std::runtime_error("Non Implemented Fast Summation Method");
    }

    std::cout << "BiotSavartExecution " << biot_savart_duration.seconds() << std::endl;

    return result;
}
