//
//  VS3DExplicit.cpp
//  MultiTracker
//
//  Created by Fang Da on 15/1/27.
//
//

#include "SimOptions.h"
#include "VS3D.h"
#include "BiotSavart.h"

namespace
{
double
angleAroundAxis(const Vec3d& v0,
                const Vec3d& v1,
                const Vec3d& a) // angle from v0 to v1 around axis a
{
    double asq = a.squaredNorm();
    assert(asq != 0);

    Vec3d u = v0 - v0.dot(a) * a / asq;
    assert(u.squaredNorm() != 0);
    u.normalize();

    Vec3d v = a.cross(u).normalized();

    return atan2(v1.dot(v), v1.dot(u));
}

}


//#define FANGS_VERSION
//#define FANGS_PATCHED

void
VS3D::step_explicit(double dt)
{
    std::cout << "Explicit time stepping" << std::endl;
    m_dbg_t1.clear();
    m_dbg_t2.clear();
    m_dbg_e1.clear();
    m_dbg_v1.clear();

    m_dbg_t1.resize(mesh().nt());
    m_dbg_t2.resize(mesh().nt());
    m_dbg_e1.resize(mesh().ne(), std::vector<double>(m_nregion, 0));
    m_dbg_v1.resize(mesh().nv(), std::vector<double>(m_nregion, 0));

    // integrate surface tension force

    std::vector<std::vector<double>> curvature(
      mesh().ne(), std::vector<double>(m_nregion, 0)); // edge-aligned curvature (signed scalar)
    std::vector<std::vector<double>> mean_curvatures(
      mesh().nv(),
      std::vector<double>(m_nregion, 0)); // vertex-aligned mean curvature (signed scalar)
    std::vector<double> avg_vertex_areas(mesh().nv(), 0.);

    for (size_t i = 0; i < mesh().nv(); i++)
    {
        std::set<int> incident_region;
        double vertex_area = 0.0;
        for (size_t j = 0; j < mesh().m_vertex_to_triangle_map[i].size(); j++)
        {
            const LosTopos::Vec3st& t = mesh().m_tris[mesh().m_vertex_to_triangle_map[i][j]];
            const LosTopos::Vec2i& l =
              mesh().get_triangle_label(mesh().m_vertex_to_triangle_map[i][j]);

            for (int r = 0; r < 2; ++r)
            {
                incident_region.insert(l[r]);
            }

            Vec3d x0 = pos(t[0]);
            Vec3d x1 = pos(t[1]);
            Vec3d x2 = pos(t[2]);
            double area = (x1 - x0).cross(x2 - x0).norm() / 2;

            vertex_area += area / 3 * 2.0;
        }

        if (incident_region.size() > 0)
        {
            avg_vertex_areas[i] = vertex_area / (double)incident_region.size();
        }
    }

    for (size_t region = 0; region < m_nregion; region++)
    {
        for (size_t i = 0; i < mesh().ne(); i++)
        {
            std::vector<size_t>
              incident_faces; // faces incident to edge i that have the label of interest (assume
                              // there are only two of them for now; this can be false only when
                              // complex collision prevents immediate T1 resolution, which is not
                              // expected to happen for bubble complexes.)
            std::set<int> incident_regions;
            for (size_t j = 0; j < mesh().m_edge_to_triangle_map[i].size(); j++)
            {
                const LosTopos::Vec2i& l =
                  mesh().get_triangle_label(mesh().m_edge_to_triangle_map[i][j]);
                if (l[0] == region || l[1] == region)
                    incident_faces.push_back(j);
                incident_regions.insert(l[0]);
                incident_regions.insert(l[1]);
            }
            if (incident_faces.size() == 0)
                continue;
            //            assert(incident_faces.size() == 2);
            if (incident_faces.size() != 2)
            {
                //                std::cout << "Warning: incident_faces.size() != 2" << std::endl;
                curvature[i][region] = 0;
                continue;
            }
            bool nonmanifold = (incident_regions.size() > 2);

            //            assert(mesh().m_edge_to_triangle_map[i].size() == 2);
            int v0 = mesh().m_edges[i][0];
            int v1 = mesh().m_edges[i][1];
            Vec3d x0 = pos(v0);
            Vec3d x1 = pos(v1);
            Vec3d et = (x1 - x0);

            int ti0 = mesh().m_edge_to_triangle_map[i][incident_faces[0]];
            int ti1 = mesh().m_edge_to_triangle_map[i][incident_faces[1]];
            LosTopos::Vec3st t0 = mesh().get_triangle(ti0);
            LosTopos::Vec3st t1 = mesh().get_triangle(ti1);

            if (mesh().get_triangle_label(ti0)[mesh().oriented(v0, v1, t0) ? 1 : 0] == region)
                std::swap(ti0, ti1),
                  std::swap(t0, t1); // the region of interest should be to the CCW direction of t0
                                     // when looking along the direciton of edge i

            assert(mesh().get_triangle_label(ti0)[mesh().oriented(v0, v1, t0) ? 0 : 1] == region);
            assert(mesh().get_triangle_label(ti1)[mesh().oriented(v0, v1, t1) ? 1 : 0] == region);

            Vec3d n0 = (pos(t0[1]) - pos(t0[0])).cross(pos(t0[2]) - pos(t0[0]));
            if (mesh().get_triangle_label(ti0)[1] == region)
                n0 = -n0; // n0 should point away from the region of interest

            Vec3d n1 = (pos(t1[1]) - pos(t1[0])).cross(pos(t1[2]) - pos(t1[0]));
            if (mesh().get_triangle_label(ti1)[1] == region)
                n1 = -n1; // n1 should point away from the region of interest

            double edgeArea = (n0.norm() + n1.norm()) / 2 / 3;
            n0.normalize();
            n1.normalize();

            //            double curvature_i = 0;
            //            if (nonmanifold)
            //            {
            //                curvature_i = n0.cross(n1).dot(et) / (et.norm() * m_delta); //
            //                integral curvature along the curve orghotonal to the edge in plane
            //            } else
            //            {
            //                curvature_i = n0.cross(n1).dot(et) / edgeArea;
            ////                curvature_i = angleAroundAxis(n0, n1, et) * et.norm() / edgeArea;
            //            }

            //            double curvature_i = n0.cross(n1).dot(et);
            //            double curvature_i = (n0 - n1).dot((n0 + n1).normalized().cross(et));
            double curvature_i = angleAroundAxis(n0, n1, et) * et.norm();

            curvature[i][region] = curvature_i;
            //            if (region == 1)
            //                m_dbg_e1[i] = curvature_i;
            m_dbg_e1[i][region] = curvature_i;
        }

        //        for (size_t i = 0; i < mesh().nv(); i++)
        //        {
        //            double mean_curvature = 0;
        //            for (size_t j = 0; j < mesh().m_vertex_to_edge_map[i].size(); j++)
        //                mean_curvature += curvature[mesh().m_vertex_to_edge_map[i][j]][region];
        //            mean_curvature /= mesh().m_vertex_to_edge_map[i].size();
        //
        //            (*m_Gamma)[i][region] += simOptions().sigma * mean_curvature * dt;
        //            m_dbg_v1[i][region] = mean_curvature;
        //        }

        for (size_t i = 0; i < mesh().nv(); i++)
        {
            double mean_curvature = 0;

            Mat3d second_fundamental_form = Mat3d::Zero();
            int counter = 0;

#ifdef FANGS_VERSION
            double vertex_area = 0;
            double triple_junction_length_sum = 0;
#else
            double vertex_area = avg_vertex_areas[i];
#endif
            for (size_t j = 0; j < mesh().m_vertex_to_edge_map[i].size(); j++)
            {
                size_t e = mesh().m_vertex_to_edge_map[i][j];
                bool incident_to_region = false;
                for (size_t k = 0; k < mesh().m_edge_to_triangle_map[e].size(); k++)
                {
                    const LosTopos::Vec2i& l =
                      mesh().get_triangle_label(mesh().m_edge_to_triangle_map[e][k]);
                    if (l[0] == region || l[1] == region)
                    {
                        incident_to_region = true;
                        break;
                    }
                }

                if (incident_to_region)
                {
                    Vec3d et = (pos(mesh().m_edges[e][1]) - pos(mesh().m_edges[e][0])).normalized();
                    second_fundamental_form += et * et.transpose() * curvature[e][region];
                    mean_curvature += curvature[e][region];
                    counter++;
#ifdef FANGS_VERSION
                    if (mesh().m_edge_to_triangle_map[e].size() > 2)
                        triple_junction_length_sum +=
                          (pos(mesh().m_edges[e][1]) - pos(mesh().m_edges[e][0])).norm();
#endif
                }
            }
#ifdef FANGS_VERSION
            for (size_t j = 0; j < mesh().m_vertex_to_triangle_map[i].size(); j++)
            {
                const LosTopos::Vec3st& t = mesh().m_tris[mesh().m_vertex_to_triangle_map[i][j]];
                const LosTopos::Vec2i& l =
                  mesh().get_triangle_label(mesh().m_vertex_to_triangle_map[i][j]);
                if (l[0] == region || l[1] == region)
                {
                    Vec3d x0 = pos(t[0]);
                    Vec3d x1 = pos(t[1]);
                    Vec3d x2 = pos(t[2]);
                    double area = (x1 - x0).cross(x2 - x0).norm() / 2;
                    vertex_area += area / 3;
                }
            }
#endif
            //            if (counter == 0)
            //                mean_curvature = 0;
            //            else
            ////                mean_curvature = second_fundamental_form.trace() / counter / 3;
            //                mean_curvature /= counter;

#ifdef FANGS_VERSION
#ifdef FANGS_PATCHED
            if (triple_junction_length_sum > 0) // this means the vertex is a triple junction
                                                // vertex; the vertex domain should then be smaller.
                vertex_area = triple_junction_length_sum * m_delta * 1.0;
#endif
#endif

            if (vertex_area == 0)
                mean_curvature = 0;
            else
                mean_curvature = mean_curvature / (vertex_area * 2);

            //            (*m_Gamma)[i][region] += simOptions().sigma * mean_curvature * dt;
            m_dbg_v1[i][region] = mean_curvature;

            mean_curvatures[i][region] = mean_curvature;
        }
    }

    // Integrate vertex mean curvatures into vertex Gammas (skipping constrained vertices)
    std::vector<bool> constrained(mesh().nv(), false);
    for (size_t i = 0; i < m_constrained_vertices.size(); i++)
        constrained[m_constrained_vertices[i]] = true;

    for (size_t i = 0; i < mesh().nv(); i++)
    {
        if (constrained[i])
            continue;

        std::set<Vec2i, Vec2iComp> incident_region_pairs;
        for (size_t j = 0; j < mesh().m_vertex_to_triangle_map[i].size(); j++)
        {
            LosTopos::Vec2i l = mesh().get_triangle_label(mesh().m_vertex_to_triangle_map[i][j]);
            incident_region_pairs.insert(l[0] < l[1] ? Vec2i(l[0], l[1]) : Vec2i(l[1], l[0]));
        }

        for (std::set<Vec2i, Vec2iComp>::iterator j = incident_region_pairs.begin();
             j != incident_region_pairs.end();
             j++)
        {
            Vec2i rp = *j;
            double mean_curvature = mean_curvatures[i][rp[0]] - mean_curvatures[i][rp[1]];
            (*m_Gamma)[i].set(rp, (*m_Gamma)[i].get(rp) + simOptions().sigma * mean_curvature * dt);
        }
    }

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
        for (const Vec2i& region_pair : getVertexIncidentRegionPairs(i))
        {
            Gamma(i).set(region_pair, Gamma(i).get(region_pair) * pow(simOptions().damping_coef, dt));

        }

    std::cout << "Explicit time stepping finished" << std::endl;
}
