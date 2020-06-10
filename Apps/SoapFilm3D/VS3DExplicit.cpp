//
//  VS3DExplicit.cpp
//  MultiTracker
//
//  Created by Fang Da on 15/1/27.
//
//

#include "BiotSavart.h"
#include "SimOptions.h"
#include "VS3D.h"

#include "fmmtl/fmmtl/util/Clock.hpp"

//#define FANGS_VERSION
//#define FANGS_PATCHED

void
VS3D::step_explicit(double dt)
{
    Clock circulation_time_stepping_duration;
    std::cout << "Explicit time stepping" << std::endl;
    m_dbg_t1.clear();
    m_dbg_t2.clear();
    m_dbg_e1.clear();
    m_dbg_v1.clear();

    m_dbg_t1.resize(mesh().nt());
    m_dbg_t2.resize(mesh().nt());
    m_dbg_v1.resize(mesh().nv());

    m_dbg_e1 = getEdgeAlignedCurvatures();
    m_dbg_v1 = getMeanCurvatures(m_dbg_e1);
    std::vector<std::map<int, double>> mean_curvatures = m_dbg_v1;



    // integrate surface tension force
    // Integrate vertex mean curvatures into vertex Gammas (skipping constrained vertices)
    std::vector<bool> constrained(mesh().nv(), false);
    for (size_t i = 0; i < m_constrained_vertices.size(); i++)
        constrained[m_constrained_vertices[i]] = true;

    for (size_t vertex_index : boost::irange(0lu, mesh().nv()))
    {
        if (constrained[vertex_index])
            continue;

        for (const Vec2i& region_pair : getVertexIncidentRegionPairs(vertex_index))
        {
            double mean_curvature = mean_curvatures[vertex_index][region_pair[0]] - mean_curvatures[vertex_index][region_pair[1]];
            Gamma(vertex_index).set(region_pair, Gamma(vertex_index).get(region_pair) + simOptions().sigma * mean_curvature * dt);
        }
    }

    std::cout << "CirculationTimeSteppingExecutionTime "
              << circulation_time_stepping_duration.seconds() << std::endl;

    // Biot-Savart integral based on vortex sheet strength gamma: Stock 2006, Eq. 2.26
    VecXd v;
    if (Options::boolValue("RK4-velocity-integration"))
    {
        // RK4 integration
        VecXd v1 = BiotSavart(*this, VecXd::Zero(mesh().nv() * 3));
        VecXd v2 = BiotSavart(*this, v1 * dt * 0.5);
        VecXd v3 = BiotSavart(*this, v2 * dt * 0.5);
        VecXd v4 = BiotSavart(*this, v3 * dt);
        v = (v1 + 2 * v2 + 2 * v3 + v4) / 6;
    }
    else
    {
        // explicit Euler
        v = BiotSavart(*this, VecXd::Zero(mesh().nv() * 3));
    }

    for (size_t i = 0; i < mesh().nv(); i++)
        m_st->pm_newpositions[i] = m_st->pm_positions[i] + vc(v.segment<3>(i * 3)) * dt;

    // damping
    for (size_t i = 0; i < mesh().nv(); i++)
    {
        for (const Vec2i& region_pair : getVertexIncidentRegionPairs(i))
        {
            Gamma(i).set(region_pair,
                         Gamma(i).get(region_pair) * pow(simOptions().damping_coef, dt));
        }
    }

    std::cout << "Explicit time stepping finished" << std::endl;
}
