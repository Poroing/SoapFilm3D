//
//  VS3D.cpp
//  MultiTracker
//
//  Created by Fang Da on 10/27/14.
//
//

#include "VS3D.h"
#include "ADT/adreal.h"
#include "ADT/advec.h"
#include "LinearizedImplicitEuler.h"
#include "LosTopos/LosTopos3D/subdivisionscheme.h"
#include "SimOptions.h"

#include "LinearBendingForce.h"
#include "SimpleGravityForce.h"
#include "SpringForce.h"
#include "VertexAreaForce.h"

#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseLU>
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>

VecXd
BiotSavart(VS3D& vs, const VecXd& dx);

bool
VS3D::isVertexConstrained(size_t vert) const
{
    return m_constrained_mapping.find(vert) != m_constrained_mapping.end();
}

VS3D::VS3D(const std::vector<LosTopos::Vec3d>& vs,
           const std::vector<LosTopos::Vec3st>& fs,
           const std::vector<LosTopos::Vec2i>& ls,
           const std::vector<size_t>& constrained_vertices,
           const std::vector<Vec3d>& constrained_positions,
           const std::vector<Vec3d>& constrained_velocities,
           const std::vector<unsigned char>& constrained_fixed)
{
    // load sim options
    m_sim_options.implicit = Options::boolValue("implicit-integration");
    m_sim_options.pbd = Options::boolValue("pbd-implicit");
    m_sim_options.smoothing_coef = Options::doubleValue("smoothing-coef");
    if (Options::strValue("smoothing-type") == "biharmonic")
    {
        m_sim_options.smoothing_type = SimOptions::SmoothingType::BIHARMONIC;
    }
    else if (Options::strValue("smoothing-type") == "umbrella")
    {
        m_sim_options.smoothing_type = SimOptions::SmoothingType::UMBRELLA;
    }
    else
    {
        m_sim_options.smoothing_type = SimOptions::SmoothingType::LAPLACIAN;
    }
    m_sim_options.damping_coef = Options::doubleValue("damping-coef");
    m_sim_options.sigma = Options::doubleValue("sigma");
    m_sim_options.gravity = Options::doubleValue("gravity");
    m_sim_options.looped = Options::boolValue("looped");
    m_sim_options.radius = Options::doubleValue("radius");
    m_sim_options.density = Options::doubleValue("density");
    m_sim_options.stretching = Options::doubleValue("stretching");
    m_sim_options.bending = Options::doubleValue("bending");

    // construct the surface tracker
    double mean_edge_len = Options::doubleValue("remeshing-resolution");
    if (mean_edge_len == 0)
    {
        for (size_t i = 0; i < fs.size(); i++)
        {
            mean_edge_len += mag(vs[fs[i][0]] - vs[fs[i][1]]);
            mean_edge_len += mag(vs[fs[i][1]] - vs[fs[i][2]]);
            mean_edge_len += mag(vs[fs[i][2]] - vs[fs[i][0]]);
        }
        mean_edge_len /= (fs.size() * 3);
    }
    double min_edge_len = mean_edge_len * 0.5;
    double max_edge_len = mean_edge_len * 1.5;

    std::cout << "mean edge length = " << mean_edge_len << " min edge length = " << min_edge_len
              << " max edge length = " << max_edge_len << std::endl;

    LosTopos::SurfTrackInitializationParameters params;
    params.m_proximity_epsilon =
      Options::doubleValue("lostopos-collision-epsilon-fraction") * mean_edge_len;
    params.m_merge_proximity_epsilon =
      Options::doubleValue("lostopos-merge-proximity-epsilon-fraction") * mean_edge_len;
    params.m_allow_vertex_movement_during_collapse = true;
    params.m_perform_smoothing = Options::boolValue("lostopos-perform-smoothing");
    params.m_min_edge_length = min_edge_len;
    params.m_max_edge_length = max_edge_len;
    params.m_max_volume_change =
      Options::doubleValue("lostopos-max-volume-change-fraction") * pow(mean_edge_len, 3);
    params.m_min_triangle_angle = Options::doubleValue("lostopos-min-triangle-angle");
    params.m_max_triangle_angle = Options::doubleValue("lostopos-max-triangle-angle");
    params.m_large_triangle_angle_to_split =
      Options::doubleValue("lostopos-large-triangle-angle-to-split");
    params.m_min_triangle_area =
      Options::doubleValue("lostopos-min-triangle-area-fraction") * pow(mean_edge_len, 2);
    params.m_verbose = false;
    params.m_allow_non_manifold = Options::boolValue("lostopos-allow-non-manifold");
    params.m_allow_topology_changes = Options::boolValue("lostopos-allow-topology-changes");
    params.m_collision_safety = true;
    params.m_remesh_boundaries = true;
    params.m_t1_transition_enabled = Options::boolValue("lostopos-t1-transition-enabled");
    params.m_pull_apart_distance =
      Options::doubleValue("lostopos-t1-pull-apart-distance-fraction") * mean_edge_len;

    params.m_velocity_field_callback = NULL;

    if (Options::boolValue("lostopos-smooth-subdivision"))
        params.m_subdivision_scheme = new LosTopos::ModifiedButterflyScheme();
    else
        params.m_subdivision_scheme = new LosTopos::MidpointScheme();

    params.m_use_curvature_when_collapsing = false;
    params.m_use_curvature_when_splitting = false;

    m_constrained_vertices = constrained_vertices;
    m_constrained_positions = constrained_positions;
    m_constrained_mapping.insert(constrained_vertices.begin(), constrained_vertices.end());

    if (m_constrained_positions.size() > 0)
    {
        m_constrained_velocities = constrained_velocities;
        if (m_constrained_velocities.size() != m_constrained_positions.size())
        {
            m_constrained_velocities.resize(m_constrained_positions.size(), Vec3d(0, 0, 0));
        }

        m_constrained_fixed = constrained_fixed;
        if (m_constrained_fixed.size() != m_constrained_positions.size())
        {
            m_constrained_fixed.resize(m_constrained_positions.size(), false);
        }

        m_constrained_mass.resize(m_constrained_positions.size() * 3);
        const size_t np = m_constrained_positions.size();
        if (m_sim_options.looped)
        {
            const size_t ne = np;
            VecXd rest_length(ne);
            for (size_t i = 0; i < ne; ++i)
            {
                int iforw = (i == np - 1) ? 0 : (i + 1);
                rest_length(i) =
                  (m_constrained_positions[iforw] - m_constrained_positions[i]).norm();
            }

            for (size_t i = 0; i < np; ++i)
            {
                int ieforw = i;
                int ieback = (i == 0) ? (ne - 1) : (i - 1);
                double len = (rest_length(ieforw) + rest_length(ieback)) * 0.5;
                double mass =
                  M_PI * m_sim_options.radius * m_sim_options.radius * m_sim_options.density * len;
                for (size_t r = 0; r < 3; ++r)
                    m_constrained_mass[i * 3 + r] = mass;
            }

            // add spring forces
            for (size_t i = 0; i < ne; ++i)
            {
                int iforw = (i == np - 1) ? 0 : (i + 1);
                m_forces.push_back(new SpringForce(
                  std::pair<int, int>(i, iforw), m_sim_options.stretching, rest_length(i)));
            }

            // add bending forces
            for (size_t i = 0; i < np; ++i)
            {
                int iforw = (i == np - 1) ? 0 : (i + 1);
                int iback = (i == 0) ? (np - 1) : (i - 1);

                int ieforw = i;
                int ieback = (i == 0) ? (ne - 1) : (i - 1);

                m_forces.push_back(new LinearBendingForce(iback,
                                                          i,
                                                          iforw,
                                                          m_sim_options.bending,
                                                          m_sim_options.bending
                                                            * m_sim_options.damping_coef * 0.0001,
                                                          Vector2s::Zero(),
                                                          rest_length(ieback),
                                                          rest_length(ieforw)));
            }
        }
        else
        {
            const size_t ne = np - 1;
            VecXd rest_length(ne);
            for (size_t i = 0; i < ne; ++i)
            {
                rest_length(i) =
                  (m_constrained_positions[i + 1] - m_constrained_positions[i]).norm();
            }
            for (size_t i = 0; i < np; ++i)
            {
                double len;
                if (i == 0)
                {
                    len = rest_length(0) * 0.5;
                }
                else if (i == np - 1)
                {
                    len = rest_length(ne - 1) * 0.5;
                }
                else
                {
                    int ieforw = i;
                    int ieback = i - 1;
                    len = (rest_length(ieforw) + rest_length(ieback)) * 0.5;
                }

                double mass =
                  M_PI * m_sim_options.radius * m_sim_options.radius * m_sim_options.density * len;
                for (size_t r = 0; r < 3; ++r)
                    m_constrained_mass[i * 3 + r] = mass;
            }

            for (size_t i = 0; i < ne; ++i)
            {
                m_forces.push_back(new SpringForce(
                  std::pair<int, int>(i, i + 1), m_sim_options.stretching, rest_length(i)));
            }

            for (size_t i = 1; i < np - 1; ++i)
            {
                int iforw = i + 1;
                int iback = i - 1;

                int ieback = i - 1;
                int ieforw = i;

                m_forces.push_back(
                  new LinearBendingForce(iback,
                                         i,
                                         iforw,
                                         m_sim_options.bending,
                                         m_sim_options.bending * m_sim_options.damping_coef,
                                         Vector2s::Zero(),
                                         rest_length(ieback),
                                         rest_length(ieforw)));
            }
        }
    }

    std::vector<LosTopos::Vec3d> masses(vs.size(), LosTopos::Vec3d(1, 1, 1));
    for (size_t i = 0; i < m_constrained_vertices.size(); i++)
        masses[m_constrained_vertices[i]] *= std::numeric_limits<double>::infinity();
    m_st = new LosTopos::SurfTrack(vs, fs, ls, masses, params);
    m_st->m_solid_vertices_callback = this;
    m_st->m_mesheventcallback = this;

    // find out the number of regions
    m_nregion = 0;
    for (size_t i = 0; i < mesh().nt(); i++)
    {
        LosTopos::Vec2i l = mesh().get_triangle_label(i);
        m_nregion = std::max(m_nregion, l[0] + 1);
        m_nregion = std::max(m_nregion, l[1] + 1);
    }

    // initialize the sheet quantities
    m_Gamma = new LosTopos::NonDestructiveTriMesh::VertexData<GammaType>(&(m_st->m_mesh));
    for (size_t i = 0; i < mesh().nv(); i++)
        (*m_Gamma)[i] = GammaType(m_nregion);

    // Biot-Savart kernel regularization parameter
    m_delta = max_edge_len * 0.5;

    m_forces.push_back(new SimpleGravityForce(Vector3s(0, 0, -m_sim_options.gravity)));
    m_forces.push_back(new VertexAreaForce(this, m_sim_options.sigma));
    // initialize constraint stepper
    m_constraint_stepper = new LinearizedImplicitEuler();
}

VS3D::VS3D(const std::vector<LosTopos::Vec3d>& vs,
           const std::vector<LosTopos::Vec3st>& fs,
           const std::vector<LosTopos::Vec2i>& ls,
           const std::vector<Vec3d>& initial_velocity_direction,
           const std::vector<double>& initial_velocity_magnitude,
           const std::vector<size_t>& constrained_vertices,
           const std::vector<Vec3d>& constrained_positions,
           const std::vector<Vec3d>& constrained_velocities,
           const std::vector<unsigned char>& constrained_fixed)
  : VS3D(vs,
         fs,
         ls,
         constrained_vertices,
         constrained_positions,
         constrained_velocities,
         constrained_fixed)
{
    projectVelocity(initial_velocity_direction, initial_velocity_magnitude);
}

VS3D::~VS3D()
{
    if (m_st)
        delete m_st;
    for (Force* f : m_forces)
        delete f;

    if (m_constraint_stepper)
        delete m_constraint_stepper;
}

std::vector<size_t>
VS3D::getNumberVerticesIncidentToRegions() const
{
    std::vector<size_t> number_vertices_incident(nregion(), 0);
    for (std::size_t vertex_index = 0; vertex_index < mesh().nv(); ++vertex_index)
    {
        for (std::size_t incident_region : getVertexIncidentRegions(vertex_index))
        {
            ++number_vertices_incident[incident_region];
        }
    }
    return number_vertices_incident;
}

double
VS3D::getBoundingBoxVolume() const
{
    double min_x = std::numeric_limits<double>::infinity();
    double min_y = std::numeric_limits<double>::infinity();
    double max_x = -min_x;
    double max_y = -min_y;
    for (std::size_t vertex_index = 0; vertex_index < mesh().nv(); ++vertex_index)
    {
        Vec3d vertex = pos(vertex_index);
        min_x = std::min(vertex.x(), min_x);
        min_y = std::min(vertex.y(), min_y);
        max_x = std::max(vertex.x(), max_x);
        max_y = std::max(vertex.y(), max_y);
    }
    return (max_x - min_x) * (max_y - min_y);
}

namespace
{
Mat3d
skewSymmetric(const Vec3d& v)
{
    Mat3d ss = Mat3d::Zero();
    ss(0, 1) = -v(2);
    ss(1, 0) = v(2);
    ss(0, 2) = v(1);
    ss(2, 0) = -v(1);
    ss(1, 2) = -v(0);
    ss(2, 1) = v(0);
    return ss;
}
}

void
VS3D::stepConstrainted(const scalar& dt)
{
    m_constraint_stepper->stepScene(*this, dt);
}

double
VS3D::step(double dt)
{
    static int counter = 0;
    counter++;

    if (counter % 2 == 0)
    {
        // mesh improvement
        improveMesh(Options::intValue("remeshing-iterations"));
    }
    else
    {
        // update gamma due to external forces (surface tension, gravity, etc), and update velocity
        // from gamma
        if (simOptions().implicit)
        {
            if (simOptions().pbd)
                step_PBD_implicit(dt);
            else
                step_implicit(dt);
        }
        else
        {
            step_explicit(dt);
        }
    }

    // shift all Gammas for each region pair to their global mean
    MatXd means(m_nregion, m_nregion);
    means.setZero();
    Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> counters(m_nregion, m_nregion);
    counters.setZero();

    for (size_t i = 0; i < mesh().nv(); i++)
    {
        std::set<Vec2i, Vec2iComp> region_pairs;
        for (size_t j = 0; j < mesh().m_vertex_to_triangle_map[i].size(); j++)
        {
            LosTopos::Vec2i l = mesh().get_triangle_label(mesh().m_vertex_to_triangle_map[i][j]);
            region_pairs.insert(l[0] < l[1] ? Vec2i(l[0], l[1]) : Vec2i(l[1], l[0]));
        }

        for (std::set<Vec2i, Vec2iComp>::iterator j = region_pairs.begin(); j != region_pairs.end();
             j++)
        {
            Vec2i rp = *j;
            means(rp[0], rp[1]) += (*m_Gamma)[i].get(rp);
            counters(rp[0], rp[1])++;
        }
    }

    for (int i = 0; i < m_nregion; i++)
        for (int j = i + 1; j < m_nregion; j++)
            if (counters(i, j) != 0)
                means(i, j) /= counters(i, j);

    for (size_t i = 0; i < mesh().nv(); i++)
    {
        std::set<Vec2i, Vec2iComp> region_pairs;
        for (size_t j = 0; j < mesh().m_vertex_to_triangle_map[i].size(); j++)
        {
            LosTopos::Vec2i l = mesh().get_triangle_label(mesh().m_vertex_to_triangle_map[i][j]);
            region_pairs.insert(l[0] < l[1] ? Vec2i(l[0], l[1]) : Vec2i(l[1], l[0]));
        }

        for (std::set<Vec2i, Vec2iComp>::iterator j = region_pairs.begin(); j != region_pairs.end();
             j++)
        {
            Vec2i rp = *j;
            (*m_Gamma)[i].set(rp, (*m_Gamma)[i].get(rp) - means(rp[0], rp[1]));
        }
    }

    // damping by smoothing
    if (true)
    {
        smoothCirculation(dt);
    }

    // before enforcing constraints, first scan through the mesh to find any solid vertices not
    // registered as constraints. they can appear due to remeshing (splitting an all-solid edge)
    std::vector<int> constrained_vertices_map(mesh().nv(), -1);
    for (size_t i = 0; i < m_constrained_vertices.size(); i++)
        constrained_vertices_map[m_constrained_vertices[i]] = i;

    /*    for (size_t i = 0; i < mesh().nv(); i++)
        {
            if (m_st->vertex_is_any_solid(i) && constrained_vertices_map[i] < 0)
            {
                m_constrained_vertices.push_back(i);
                m_constrained_positions.push_back(pos(i));
            }
        }*/

    // contruct the open boudnary extra faces
    m_obefc.clear();
    m_obefe.clear();
    m_obefn.clear();
    std::vector<size_t> ob;   // open boundary edges
    std::set<size_t> obv_set; // open boundary vertices
    for (size_t i = 0; i < mesh().ne(); i++)
    {
        if (surfTrack()->edge_is_all_solid(i)) // this is a constrained edge
        {
            // assert(mesh().m_edge_to_triangle_map[i].size() == 2);
            if (mesh().m_edge_to_triangle_map[i].size() == 2)
            {

                size_t f0 = mesh().m_edge_to_triangle_map[i][0];
                size_t f1 = mesh().m_edge_to_triangle_map[i][1];
                bool s0 = surfTrack()->triangle_is_all_solid(f0);
                bool s1 = surfTrack()->triangle_is_all_solid(f1);

                if ((s0 && !s1) || (s1 && !s0)) // this is an open boundary edge
                {
                    size_t f = (s0 ? f1 : f0);

                    size_t v0 = mesh().m_edges[i][0];
                    size_t v1 = mesh().m_edges[i][1];
                    size_t v2 = mesh().get_third_vertex(i, f);

                    Vec3d x0 = pos(v0);
                    Vec3d x1 = pos(v1);
                    Vec3d x2 = pos(v2);

                    m_obefc.push_back((x0 + x1) / 2);
                    m_obefe.push_back((x1 - x0).normalized());
                    m_obefn.push_back((x1 - x0).cross(x2 - x0).normalized());
                    ob.push_back(i);
                    obv_set.insert(v0);
                    obv_set.insert(v1);
                }
            }
            else if (mesh().m_edge_to_triangle_map[i].size() == 1)
            {
                size_t f = mesh().m_edge_to_triangle_map[i][0];

                size_t v0 = mesh().m_edges[i][0];
                size_t v1 = mesh().m_edges[i][1];
                size_t v2 = mesh().get_third_vertex(i, f);

                Vec3d x0 = pos(v0);
                Vec3d x1 = pos(v1);
                Vec3d x2 = pos(v2);

                m_obefc.push_back((x0 + x1) / 2);
                m_obefe.push_back((x1 - x0).normalized());
                m_obefn.push_back((x1 - x0).cross(x2 - x0).normalized());
                ob.push_back(i);
                obv_set.insert(v0);
                obv_set.insert(v1);
            }
        }
    }

    std::vector<size_t> obv;
    obv.assign(obv_set.begin(), obv_set.end());
    assert(obv.size() == ob.size()); // assume the open boundary has simple topology
    size_t nob = ob.size();

    if (nob > 0)
    {

        std::vector<Vec3d> obvn(nob); // open boundary vertex normals
        for (size_t i = 0; i < nob; i++)
        {
            Vec3d n(0, 0, 0);
            for (size_t j = 0; j < nob; j++)
                if (mesh().m_edges[ob[j]][0] == obv[i] || mesh().m_edges[ob[j]][1] == obv[i])
                    n += m_obefn[j];
            obvn[i] = n.normalized();
        }

        // solve for the obef vorticities such that the open boundary does not move
        MatXd obA = MatXd::Zero(nob, nob);
        for (size_t i = 0; i < nob; i++)
        {
            Vec3d x = pos(obv[i]);

            for (size_t j = 0; j < nob; j++)
            {
                Vec3d xp = m_obefc[j];

                Vec3d dx = x - xp;
                double dxn = sqrt(dx.squaredNorm() + m_delta * m_delta);

                obA(i, j) = -(skewSymmetric(dx) * m_obefe[j]).dot(obvn[i]) / (dxn * dxn * dxn);
            }
        }

        obA *= dt / (4 * M_PI);

        VecXd obrhs = VecXd::Zero(nob);
        for (size_t i = 0; i < nob; i++)
            obrhs[i] = (m_constrained_positions[constrained_vertices_map[obv[i]]]
                        - vc(surfTrack()->pm_newpositions[obv[i]]))
                         .dot(obvn[i]);

        //    VecXd obefv = obA.partialPivLu().solve(obrhs);
        double oblambda = 0.1;
        VecXd obefv = (obA.transpose() * obA + oblambda * oblambda * MatXd::Identity(nob, nob))
                        .partialPivLu()
                        .solve(obA.transpose() * obrhs); // regularized solve to avoid blowing up in
                                                         // presence of near-dependent constraints

        m_obefv.resize(nob, 0);
        for (size_t i = 0; i < nob; i++)
            m_obefv[i] = obefv[i];
    }
    // following the open boundary solve, recompute the velocities
    VecXd newv = BiotSavart(*this, VecXd::Zero(mesh().nv() * 3));
    for (size_t i = 0; i < mesh().nv(); i++)
        m_st->pm_newpositions[i] = m_st->pm_positions[i] + vc(newv.segment<3>(i * 3)) * dt;

    // project to remove motion on the constrained vertices

    // assume that constrained vertices can only be manifold for now
    // first of all, exclude from the solve those constrained vertices being surrounded by all
    // constrained vertices (i.e. with no incident face that is not fully constrained). These
    // vertices don't contribute vorticity because
    //  their incident faces don't contribute vorticity, and they don't desire displacement
    //  correction because they just passively move to wherever they are prescribed to go. So they
    //  don't appear in either columns or rows in the solve.
    std::vector<int> constrained_vertex_nonconstrained_neighbors(
      m_constrained_vertices.size(),
      0); // how many non-fully-constrained incident faces each constrained vertex has
    for (size_t i = 0; i < m_constrained_vertices.size(); i++)
        for (size_t j = 0; j < mesh().m_vertex_to_triangle_map[m_constrained_vertices[i]].size();
             j++)
        {
            LosTopos::Vec3st t =
              mesh().get_triangle(mesh().m_vertex_to_triangle_map[m_constrained_vertices[i]][j]);
            if (!(m_st->vertex_is_any_solid(t[0]) && m_st->vertex_is_any_solid(t[1])
                  && m_st->vertex_is_any_solid(t[2])))
                constrained_vertex_nonconstrained_neighbors[i]++;
        }

    std::vector<size_t>
      relevant_constrained_vertices; // constrained vertices that are not fully constrained
                                     // (circulations on those vertices don't matter) and are not on
                                     // open boundary (the constraint solve there is different and
                                     // they have already been handled by the OB solve above)
    for (size_t i = 0; i < m_constrained_vertices.size(); i++)
    {
        if (constrained_vertex_nonconstrained_neighbors[i] == 0)
            continue; // ignore constrained vertices who don't have any non-fully-constrianed
                      // incident faces

        std::set<Vec2i, Vec2iComp> rps;
        for (size_t j = 0; j < mesh().m_vertex_to_triangle_map[m_constrained_vertices[i]].size();
             j++)
        {
            LosTopos::Vec3st t =
              mesh().get_triangle(mesh().m_vertex_to_triangle_map[m_constrained_vertices[i]][j]);
            if (m_st->vertex_is_any_solid(t[0]) && m_st->vertex_is_any_solid(t[1])
                && m_st->vertex_is_any_solid(t[2]))
                continue;
            LosTopos::Vec2i l = mesh().get_triangle_label(
              mesh().m_vertex_to_triangle_map[m_constrained_vertices[i]][j]);
            Vec2i rp = (l[0] < l[1] ? Vec2i(l[0], l[1]) : Vec2i(l[1], l[0]));
            rps.insert(rp);
        }
        if (rps.size() > 1)
            continue; // ignore constrained vertices who are incident to more than one region pair
                      // (through non-fully-constrained incident faces), because the code below
                      // can't handle them

        bool ob = false;
        for (size_t j = 0; j < nob; j++)
            if (obv[j] == m_constrained_vertices[i])
                ob = true;
        if (ob)
            continue; // ignore open boundary vertices

        relevant_constrained_vertices.push_back(i);
    }

    size_t nc = relevant_constrained_vertices.size();

#warning There should be a way to factorize this with the BiotSavart function.
    if (nc > 0)
    {
        std::vector<Vec3d> constrained_vertex_normal(nc); // surface normals
        std::vector<size_t> relevant_constrained_vertices_global_indices(nc);
        std::vector<double> target_normal_velocity;
        for (size_t i = 0; i < nc; i++)
        {
            size_t cv = m_constrained_vertices[relevant_constrained_vertices[i]];
            relevant_constrained_vertices_global_indices[i] = cv;
            constrained_vertex_normal[i] = getManifoldVertexNormal(cv);
            target_normal_velocity[i] =
              (m_constrained_positions[relevant_constrained_vertices[i]] - pos(i))
                .dot(constrained_vertex_normal[i])
              / dt;
        }
        projectVelocity(relevant_constrained_vertices_global_indices,
                        constrained_vertex_normal,
                        target_normal_velocity);
    }

#warning There might be a way to reuse the computation of the previous BiotSavart
    // recompute the velocity after the constraint projection
    newv = BiotSavart(*this, VecXd::Zero(mesh().nv() * 3));
    for (size_t i = 0; i < mesh().nv(); i++)
        m_st->pm_newpositions[i] = m_st->pm_positions[i] + vc(newv.segment<3>(i * 3)) * dt;

    // enforce the constraints exactly, potentially sliding them tangentially. TODO: resample
    // Gammas?
    for (size_t ii = 0; ii < m_constrained_vertices.size(); ii++)
    {
        size_t i = m_constrained_vertices[ii];
        m_st->pm_newpositions[i] = vc(m_constrained_positions[ii]);
    }

    // move the mesh
    if (counter % 2 != 0)
    {
        double actual_dt;
        m_st->integrate(dt, actual_dt);
        if (actual_dt != dt)
            std::cout
              << "Warning: SurfTrack::integrate() failed to step the full length of the time step!"
              << std::endl;
    }

    //    update_dbg_quantities();

    return (counter % 2 == 0 ? 0 : dt);
}

void
VS3D::smoothCirculation(double dt)
{
    if (simOptions().smoothing_type == SimOptions::SmoothingType::BIHARMONIC)
    {
        biharmonicSmoothing(dt);
    }
    else if (simOptions().smoothing_type == SimOptions::SmoothingType::UMBRELLA)
    {
        umbrellaSmoothing(dt);
    }
    else
    {
        laplacianSmoothing(dt);
    }
}

void
VS3D::biharmonicSmoothing(double dt)
{
    std::vector<Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic>> incident_region_pairs =
      getIncidentRegions();

    MatXd igl_ready_positions = getIglReadyPositions();
    MatXi igl_ready_triangles = getIglReadyTriangles();

    SparseMatd mass_matrix;
    SparseMatd laplace_operator;
    igl::cotmatrix(igl_ready_positions, igl_ready_triangles, laplace_operator);
    igl::massmatrix(
      igl_ready_positions, igl_ready_triangles, igl::MASSMATRIX_TYPE_VORONOI, mass_matrix);

    Eigen::SimplicialLDLT<SparseMatd> mass_matrix_solver(mass_matrix);
    SparseMatd inverse_mass_matrix_laplace_operator = mass_matrix_solver.solve(laplace_operator);
    SparseMatd biharmonic_energy_hessian =
      laplace_operator.transpose() * inverse_mass_matrix_laplace_operator;

    Eigen::SimplicialLLT<SparseMatd> smoothing_solver(
      (1 - simOptions().smoothing_coef) * dt * mass_matrix
      + simOptions().smoothing_coef * dt * inverse_mass_matrix_laplace_operator);

    VecXd new_gammas;
    for (int j = 0; j < m_nregion; j++)
    {
        for (int k = j + 1; k < m_nregion; k++)
        {
            new_gammas = smoothing_solver.solve((1 - simOptions().smoothing_coef) * dt * mass_matrix
                                                * getGammas(j, k));
            for (std::size_t vertex_index : boost::irange(0lu, mesh().nv()))
            {
                if (!incident_region_pairs[vertex_index](j, k))
                {
                    new_gammas[vertex_index] = 0;
                }
            }
            setGammas(j, k, new_gammas);
        }
    }
}

void
VS3D::umbrellaSmoothing(double dt)
{
    std::vector<Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic>> incident_region_pairs =
      getIncidentRegions();

    MatXd igl_ready_positions = getIglReadyPositions();
    MatXi igl_ready_triangles = getIglReadyTriangles();

    VecXd umbrella_laplacian;
    VecXd gammas;
    for (int j = 0; j < m_nregion; j++)
    {
        for (int k = j + 1; k < m_nregion; k++)
        {
            umbrella_laplacian = VecXd::Zero(mesh().nv());
            gammas = getGammas(j, k);
            for (std::size_t vertex_index : boost::irange(0lu, mesh().nv()))
            {
                if (!incident_region_pairs[vertex_index](j, k))
                {
                    continue;
                }

                for (std::size_t adjacent_vertex_index : getVertexAdjacentVertices(vertex_index))
                {
                    if (!incident_region_pairs[vertex_index](j, k))
                    {
                        continue;
                    }

                    umbrella_laplacian[vertex_index] +=
                      gammas[adjacent_vertex_index] - gammas[vertex_index];
                }

                umbrella_laplacian[vertex_index] /=
                  getVertexAreaIncidentToRegionPair(vertex_index, Vec2i(j, k));
            }
            setGammas(j, k, gammas + simOptions().smoothing_coef * dt * umbrella_laplacian);
        }
    }
}

void
VS3D::laplacianSmoothing(double dt)
{
    std::vector<Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic>> incident_region_pairs =
      getIncidentRegions();

    MatXd igl_ready_positions = getIglReadyPositions();
    MatXi igl_ready_triangles = getIglReadyTriangles();

    SparseMatd mass_matrix;
    SparseMatd laplace_operator;

    igl::massmatrix(
      igl_ready_positions, igl_ready_triangles, igl::MASSMATRIX_TYPE_VORONOI, mass_matrix);

    VecXd new_gammas;
    for (int j = 0; j < m_nregion; j++)
    {
        for (int k = j + 1; k < m_nregion; k++)
        {
            igl_ready_triangles = getIglReadyTrianglesIncidentToRegionPair(Vec2i(j, k));
            igl::cotmatrix(igl_ready_positions, igl_ready_triangles, laplace_operator);
            Eigen::SimplicialLLT<SparseMatd> smoothing_solver(
              mass_matrix - simOptions().smoothing_coef * dt * laplace_operator);

            new_gammas = smoothing_solver.solve(mass_matrix * getGammas(j, k));
            setGammas(j, k, new_gammas);
        }
    }
}

std::vector<Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic>>
VS3D::getIncidentRegions() const
{
    std::vector<Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic>> incident_region_pairs(
      mesh().nv());
    for (size_t i = 0; i < mesh().nv(); i++)
        incident_region_pairs[i].setZero(m_nregion, m_nregion);

    for (size_t i = 0; i < mesh().nt(); i++)
    {
        LosTopos::Vec3st t = mesh().get_triangle(i);
        LosTopos::Vec2i l = mesh().get_triangle_label(i);
        incident_region_pairs[t[0]](l[0], l[1]) = incident_region_pairs[t[0]](l[1], l[0]) = true;
        incident_region_pairs[t[1]](l[0], l[1]) = incident_region_pairs[t[1]](l[1], l[0]) = true;
        incident_region_pairs[t[2]](l[0], l[1]) = incident_region_pairs[t[2]](l[1], l[0]) = true;
    }
    return incident_region_pairs;
}

void
VS3D::getRegionPairIncidentTriangles(const Vec2i& region_pair,
                                     std::vector<size_t>& triangle_indices) const
{
    auto region_pair_incident_triangles = getRegionPairIncidentTriangles(region_pair);
    triangle_indices.assign(region_pair_incident_triangles.begin(),
                            region_pair_incident_triangles.end());
}

bool
VS3D::isTriangleIncidentToRegionPair(size_t triangle_index, const Vec2i& region_pair) const
{
    LosTopos::Vec2i triangle_label = mesh().get_triangle_label(triangle_index);
    return (triangle_label[0] == region_pair[0] && triangle_label[1] == region_pair[1])
           || (region_pair[0] == triangle_label[1] && region_pair[1] == triangle_label[0]);
}

size_t
VS3D::getEdgeOtherVertex(size_t edge_index, size_t vertex_index) const
{
    LosTopos::Vec2st e = mesh().m_edges[edge_index];
    return (e[0] == vertex_index ? e[1] : e[0]);
}

size_t
VS3D::getVertexDegree(size_t vertex_index) const
{
    return mesh().m_vertex_to_edge_map[vertex_index].size();
}

double
VS3D::getVertexAreaIncidentToRegionPair(size_t vertex_index, const Vec2i& region_pair) const
{
    double area = 0;
    for (size_t triangle_index : mesh().m_vertex_to_triangle_map[vertex_index])
    {
        if (!isTriangleIncidentToRegionPair(triangle_index, region_pair))
        {
            continue;
        }

        area += m_st->get_triangle_area(triangle_index) / 3;
    }
    return area;
}

MatXd
VS3D::getIglReadyPositions() const
{
    MatXd igl_ready_positions(mesh().nv(), 3u);
    for (size_t vertex_index = 0; vertex_index < mesh().nv(); ++vertex_index)
    {
        igl_ready_positions.row(vertex_index) = pos(vertex_index);
    }
    return igl_ready_positions;
}

MatXi
VS3D::getIglReadyTrianglesIncidentToRegionPair(const Vec2i& region_pair) const
{
    std::vector<size_t> triangle_incident_to_region_pair_indices;
    getRegionPairIncidentTriangles(region_pair, triangle_incident_to_region_pair_indices);
    return getIglReadyTriangles(triangle_incident_to_region_pair_indices);
}

MatXi
VS3D::getIglReadyTriangles() const
{
    MatXi igl_ready_triangles(mesh().nt(), 3u);
    for (size_t triangle_index = 0; triangle_index < mesh().nt(); ++triangle_index)
    {
        igl_ready_triangles.row(triangle_index) = vc(mesh().get_triangle(triangle_index));
    }
    return igl_ready_triangles;
}

MatXi
VS3D::getIglReadyTriangles(const std::vector<size_t>& triangle_indices) const
{
    MatXi igl_ready_triangles(triangle_indices.size(), 3u);
    for (size_t i : boost::irange(0lu, triangle_indices.size()))
    {
        igl_ready_triangles.row(i) = vc(mesh().get_triangle(triangle_indices[i]));
    }
    return igl_ready_triangles;
}

VecXd
VS3D::getGammas(std::vector<std::pair<size_t, Vec2i>>& gamma_to_vertex_and_region_pair,
                std::vector<Eigen::MatrixXi>& vertex_and_region_pair_to_gamma) const
{
    std::vector<double> gammas;
    gamma_to_vertex_and_region_pair.clear();
    vertex_and_region_pair_to_gamma.resize(mesh().nv(),
                                           -Eigen::MatrixXi::Ones(nregion(), nregion()));

    size_t gamma_index = 0;
    for (size_t vertex_index : boost::irange(0lu, mesh().nv()))
    {
        for (Vec2i region_pair : getVertexIncidentRegionPairs(vertex_index))
        {
            gamma_to_vertex_and_region_pair.push_back(std::make_pair(vertex_index, region_pair));
            vertex_and_region_pair_to_gamma[vertex_index](region_pair.x(), region_pair.y()) =
              gamma_index;
            vertex_and_region_pair_to_gamma[vertex_index](region_pair.y(), region_pair.x()) =
              gamma_index;
            gammas.push_back(Gamma(vertex_index).get(region_pair));
            ++gamma_index;
        }
    }
    return VecXd::Map(gammas.data(), gammas.size());
}

VecXd
VS3D::getGammas(size_t region_a_index, size_t region_b_index) const
{
    return getGammas(Vec2i(region_a_index, region_b_index));
}

VecXd
VS3D::getGammas(const Vec2i& region_pair) const
{
    VecXd gammas(mesh().nv());
    for (size_t vertex_index = 0; vertex_index < mesh().nv(); ++vertex_index)
    {
        gammas[vertex_index] = Gamma(vertex_index).get(region_pair);
    }
    return gammas;
}

void
VS3D::setGammas(const std::vector<std::pair<size_t, Vec2i>>& gamma_to_vertex_and_region_pair,
                const VecXd& gammas)
{
    for (size_t gamma_index : boost::irange(gammas.size()))
    {
        auto [vertex_index, region_pair] = gamma_to_vertex_and_region_pair[gamma_index];
        Gamma(vertex_index).set(region_pair, gammas[gamma_index]);
    }
}

void
VS3D::setGammas(size_t region_a_index, size_t region_b_index, const VecXd& gammas)
{
    return setGammas(Vec2i(region_a_index, region_b_index), gammas);
}

void
VS3D::setGammas(const Vec2i& region_pair, const VecXd& gammas)
{
    for (size_t vertex_index = 0; vertex_index < mesh().nv(); ++vertex_index)
    {
        Gamma(vertex_index).set(region_pair, gammas[vertex_index]);
    }
}

void
VS3D::improveMesh(size_t number_iteration)
{
    for (int i = 0; i < number_iteration; i++)
    {
        m_st->topology_changes();
        m_st->improve_mesh();
    }

    // defrag the mesh in the end, to ensure the next step starts with a clean mesh
    m_st->defrag_mesh_from_scratch(m_constrained_vertices);
    for (size_t i = 0; i < m_constrained_vertices.size(); i++)
        assert(m_constrained_vertices[i] < mesh().nv());

    std::cout << "NumberVerticesIncidentToRegions";
    for (std::size_t number_vertices : getNumberVerticesIncidentToRegions())
    {
        std::cout << " " << number_vertices;
    }
    std::cout << std::endl;
}

void
VS3D::projectVelocity(const std::vector<Vec3d>& direction,
                      const std::vector<double>& velocity_along_direction)
{
    auto vertices_indices_range = boost::irange(0lu, mesh().nv());
    std::vector<size_t> vertices_indices(vertices_indices_range.begin(),
                                         vertices_indices_range.end());
    projectVelocity(vertices_indices, direction, velocity_along_direction);
}

void
VS3D::projectVelocity(const std::vector<size_t>& vertices_indices,
                      const std::vector<Vec3d>& direction,
                      const std::vector<double>& velocity_along_direction)
{
    assert(vertices_indices.size() == direction.size());
    assert(vertices_indices.size() == velocity_along_direction.size());
    MatXd A = getCirculationToProjectedVelocityMatrix(vertices_indices, direction);

    VecXd rhs = VecXd::Zero(vertices_indices.size());
    for (size_t i = 0; i < velocity_along_direction.size(); ++i)
    {
        rhs(i) = velocity_along_direction[i];
    }

    double lambda = (m_st->m_min_edge_length + m_st->m_max_edge_length) / 2 * 0.1;
    VecXd result =
      (A.transpose() * A
       + lambda * lambda * MatXd::Identity(vertices_indices.size(), vertices_indices.size()))
        .partialPivLu()
        .solve(A.transpose() * rhs); // regularized solve to avoid blowing up in
                                     // presence of near-dependent constraints

    for (size_t i = 0; i < vertices_indices.size(); ++i)
    {
        size_t vertex_index = vertices_indices[i];
        Vec2i incident_regions = getManifoldVertexRegionPair(vertex_index);
        (*m_Gamma)[i].set(incident_regions, result[i]);
    }
}

void
VS3D::projectAirVelocity(const std::vector<Vec3d>& velocity, const std::vector<Vec3d>& positions)
{
    assert(velocity.size() == positions.size());

    std::vector<Eigen::MatrixXi> vertex_and_region_pair_to_gamma;
    std::vector<std::pair<size_t, Vec2i>> gamma_to_vertex_and_region_pair;
    VecXd gammas = getGammas(gamma_to_vertex_and_region_pair, vertex_and_region_pair_to_gamma);
    assert(velocity.size() * 3 <= gammas.size());

    size_t system_size = gammas.size() + velocity.size() * 3;

    VecXd right_hand_size(system_size);
    std::vector<Eigen::Triplet<double>> triplets;
    for (size_t gamma_index : boost::irange((VecXd::Index)0, gammas.size()))
    {
        auto [vertex_index, region_pair] = gamma_to_vertex_and_region_pair[gamma_index];
        double vertex_area_incident_to_region_pair =
          getVertexAreaIncidentToRegionPair(vertex_index, region_pair);
        triplets.emplace_back(gamma_index, gamma_index, vertex_area_incident_to_region_pair);
        right_hand_size[gamma_index] = vertex_area_incident_to_region_pair * gammas[gamma_index];
    }

    MatXd biot_savart_matrix = getSheetStrengthToVelocityMatrix(positions);
    SparseMatd circulation_to_sheet_strength_matrix = getCirculationToSheetStrengthMatrix(
      gamma_to_vertex_and_region_pair, vertex_and_region_pair_to_gamma);
    MatXd circulation_to_velocity_matrix =
      biot_savart_matrix * circulation_to_sheet_strength_matrix;
    for (size_t i : boost::irange((MatXd::Index)0, circulation_to_velocity_matrix.rows()))
    {
        for (size_t j : boost::irange((MatXd::Index)0, circulation_to_velocity_matrix.cols()))
        {
            triplets.emplace_back(gammas.size() + i, j, circulation_to_velocity_matrix(i, j));
            triplets.emplace_back(j, gammas.size() + i, circulation_to_velocity_matrix(i, j));
        }
    }

    for (size_t velocity_index : boost::irange(0lu, velocity.size()))
    {
        right_hand_size.segment<3>(gammas.size() + velocity_index * 3) = velocity[velocity_index];
    }


    std::cout << "Starting Solve" << std::endl;
    SparseMatd system_matrix(system_size, system_size);
    system_matrix.setFromTriplets(triplets.begin(), triplets.end());
    system_matrix.makeCompressed();
    Eigen::SparseLU<SparseMatd> solver;
    solver.analyzePattern(system_matrix);
    solver.factorize(system_matrix);
    VecXd result = solver.solve(right_hand_size);
    gammas = result.segment(0, gammas.size());
    std::cout << "End Solve" << std::endl;

    setGammas(gamma_to_vertex_and_region_pair, gammas);
}

MatXd
VS3D::getSheetStrengthToVelocityMatrix(const std::vector<Vec3d>& positions) const
{
    MatXd matrix(positions.size() * 3, mesh().nt() * 3);
    for (size_t velocity_index : boost::irange(0lu, positions.size()))
    {
        for (size_t triangle_index : boost::irange(0lu, mesh().nt()))
        {
            Vec3d triangle_centre = getTriangleCenter(triangle_index);
            Vec3d grad_G = getGreensFunctionGradient(positions[velocity_index]
                                                     - getTriangleCenter(triangle_index));
            matrix.block<3, 3>(velocity_index * 3, triangle_index * 3) = - skewSymmetric(grad_G);
        }
    }
    return matrix;
}

SparseMatd
VS3D::getCirculationToSheetStrengthMatrix(
  const std::vector<std::pair<size_t, Vec2i>>& gamma_to_vertex_and_region_pair,
  const std::vector<Eigen::MatrixXi>& vertex_and_region_pair_to_gamma) const
{
    std::vector<Eigen::Triplet<double>> triplets;
    for (size_t triangle_index : boost::irange(0lu, mesh().nt()))
    {
        LosTopos::Vec3st triangle = mesh().get_triangle(triangle_index);
        LosTopos::Vec2i region_pair = mesh().get_triangle_label(triangle_index);

        double sign = 1.;
        if (region_pair[0] > region_pair[1])
        {
            sign = -1.;
        }

        size_t gamma_0 =
          vertex_and_region_pair_to_gamma[triangle[0]](region_pair[0], region_pair[1]);
        size_t gamma_1 =
          vertex_and_region_pair_to_gamma[triangle[1]](region_pair[0], region_pair[1]);
        size_t gamma_2 =
          vertex_and_region_pair_to_gamma[triangle[2]](region_pair[0], region_pair[1]);
        Vec3d edge_01 = pos(triangle[1]) - pos(triangle[0]);
        Vec3d edge_12 = pos(triangle[2]) - pos(triangle[1]);
        Vec3d edge_20 = pos(triangle[0]) - pos(triangle[2]);

        triplets.emplace_back(triangle_index * 3 + 0, gamma_0, -sign * (edge_12[0]));
        triplets.emplace_back(triangle_index * 3 + 1, gamma_0, -sign * (edge_12[1]));
        triplets.emplace_back(triangle_index * 3 + 2, gamma_0, -sign * (edge_12[2]));

        triplets.emplace_back(triangle_index * 3 + 0, gamma_1, -sign * (edge_20[0]));
        triplets.emplace_back(triangle_index * 3 + 1, gamma_1, -sign * (edge_20[1]));
        triplets.emplace_back(triangle_index * 3 + 2, gamma_1, -sign * (edge_20[2]));

        triplets.emplace_back(triangle_index * 3 + 0, gamma_2, -sign * (edge_01[0]));
        triplets.emplace_back(triangle_index * 3 + 1, gamma_2, -sign * (edge_01[1]));
        triplets.emplace_back(triangle_index * 3 + 2, gamma_2, -sign * (edge_01[2]));
    }

    SparseMatd matrix(3 * mesh().nt(), gamma_to_vertex_and_region_pair.size());
    matrix.setFromTriplets(triplets.begin(), triplets.end());
    return matrix;
}

MatXd
VS3D::getCirculationToProjectedVelocityMatrix(
  const std::vector<size_t>& vertices_indices,
  const std::vector<Vec3d>& projected_velocity_direction) const
{
    assert(projected_velocity_direction.size() == vertices_indices.size());
    MatXd projection = MatXd::Zero(vertices_indices.size(), vertices_indices.size() * 3u);
    for (size_t i = 0; i < vertices_indices.size(); ++i)
    {
        projection.row(i).segment<3u>(i * 3u) = projected_velocity_direction[i];
    }
    return projection * getCirculationToVelocityMatrix(vertices_indices);
}

MatXd
VS3D::getCirculationToVelocityMatrix(const std::vector<size_t>& vertices_indices) const
{
    MatXd A = MatXd::Zero(vertices_indices.size() * 3u, vertices_indices.size());

    for (size_t i = 0; i < vertices_indices.size(); ++i)
    {
        // Biot-Savart for constrained vertex i
        size_t vertex_index = vertices_indices[i];
        Vec3d x = pos(vertex_index);

        for (size_t j = 0; j < vertices_indices.size(); j++)
        {
            size_t other_vertex_index = vertices_indices[j];
            for (size_t triangle_index : mesh().m_vertex_to_triangle_map[j])
            {
                if (surfTrack()->triangle_is_all_solid(triangle_index))
                    continue; // all-solid faces don't contribute vorticity.

                LosTopos::Vec2i l = mesh().get_triangle_label(triangle_index);
                Vec3d e_opposite = getVertexOppositeEdgeInTriangle(j, triangle_index);
                Vec3d xp = getTriangleCenter(triangle_index);

                Vec3d dx = x - xp;
                double dxn = sqrt(dx.squaredNorm() + m_delta * m_delta);
                Vec3d triangle_contribution =
                  (l[0] < l[1] ? 1 : -1) * (skewSymmetric(dx) * e_opposite / (dxn * dxn * dxn));

                A.col(j).segment<3u>(3 * i) += triangle_contribution;
            }
        }
    }
    A *= 1 / (4 * M_PI);
    return A;
}

Vec3d
VS3D::getTriangleCenter(size_t triangle_index) const
{
    LosTopos::Vec3st t = mesh().get_triangle(triangle_index);
    return (pos(t[0]) + pos(t[1]) + pos(t[2])) / 3;
}

Vec3d
VS3D::getVertexOppositeEdgeInTriangle(size_t vertex_index, size_t triangle_index) const
{
    LosTopos::Vec3st triangle = mesh().get_triangle(triangle_index);
    return pos(triangle[(vertex_index + 2) % 3]) - pos(triangle[(vertex_index + 1) % 3]);
}

Vec3d
VS3D::getManifoldVertexNormal(size_t vertex_index) const
{
    // Since we need to ignore solid triangles we can't use the
    // get_vertex_normal_angleweighted_by_label method from the surface tracking library.
    Vec3d normal(0, 0, 0);
    for (size_t triangle_index : mesh().m_vertex_to_triangle_map[vertex_index])
    {
        if (m_st->triangle_is_all_solid(triangle_index))
            continue;
        LosTopos::Vec2i region_pair = mesh().get_triangle_label(triangle_index);
        normal +=
          surfTrack()->get_triangle_area(triangle_index)
          * getTriangleNormalTowardRegion(triangle_index, std::max(region_pair[0], region_pair[1]));
    }
    assert(normal.norm() != 0);
    return normal.normalized();
}

Vec2i
VS3D::getManifoldVertexRegionPair(size_t vertex_index) const
{
    std::vector<int> incident_regions = getVertexIncidentRegions(vertex_index);
    // this vertex can't be incident to no unconstrained triangle and is manifold
    assert(incident_regions.size() == 2);
    std::sort(incident_regions.begin(), incident_regions.end());
    return Vec2i(incident_regions[0], incident_regions[1]);
}

std::vector<int>
VS3D::getVertexIncidentRegions(size_t vertex_index) const
{
    std::unordered_set<int> result;
    LosTopos::Vec2i label;
    for (std::size_t triangle_index : mesh().m_vertex_to_triangle_map[vertex_index])
    {
        if (m_st->triangle_is_all_solid(triangle_index))
        {
            continue;
        }

        label = mesh().get_triangle_label(triangle_index);
        result.insert(label[0]);
        result.insert(label[1]);
    }
    return std::vector<int>(result.begin(), result.end());
}

std::vector<Vec2i>
VS3D::getVertexIncidentRegionPairs(size_t vertex_index) const
{
    std::vector<Vec2i> incident_region_pairs;
    std::vector<int> incident_regions = getVertexIncidentRegions(vertex_index);
    for (size_t region_a_index : boost::irange(0lu, incident_regions.size()))
    {
        for (size_t region_b_index : boost::irange(region_a_index + 1, incident_regions.size()))
        {
            size_t region_a = region_a_index;
            size_t region_b = region_b_index;
            incident_region_pairs.emplace_back(std::max(region_a, region_b),
                                               std::min(region_a, region_b));
        }
    }
    return incident_region_pairs;
}

Vec3d
VS3D::getTriangleSheetStrength(size_t triangle_index) const
{
    LosTopos::Vec3st triangle = mesh().get_triangle(triangle_index);
    LosTopos::Vec2i incident_regions = mesh().get_triangle_label(triangle_index);
    Vec3d e01 = pos(triangle[1]) - pos(triangle[0]);
    Vec3d e12 = pos(triangle[2]) - pos(triangle[1]);
    Vec3d e20 = pos(triangle[0]) - pos(triangle[2]);

    return -(e01 * (*m_Gamma)[triangle[2]].get(incident_regions)
             + e12 * (*m_Gamma)[triangle[0]].get(incident_regions)
             + e20 * (*m_Gamma)[triangle[1]].get(incident_regions));
}

Vec3d
VS3D::getGreensFunctionGradient(const Vec3d& argument) const
{
    double argument_norm = std::sqrt(argument.squaredNorm() + delta() * delta());
    return -argument / (4 * M_PI * argument_norm * argument_norm * argument_norm);
}

void
VS3D::update_dbg_quantities()
{
    // compute the vertex velocities for rendering
    if (false)
    {
        m_dbg_v2.clear();
        m_dbg_v2.resize(mesh().nv());
        for (size_t i = 0; i < mesh().nv(); i++)
        {
            Vec3d v(0, 0, 0);
            Vec3d x = pos(i);

            for (size_t j = 0; j < mesh().nt(); j++)
            {
                LosTopos::Vec3st t = mesh().get_triangle(j);
                LosTopos::Vec2i l = mesh().get_triangle_label(j);
                Vec3d xp = (pos(t[0]) + pos(t[1]) + pos(t[2])) / 3;

                Vec3d e01 = pos(t[1]) - pos(t[0]);
                Vec3d e12 = pos(t[2]) - pos(t[1]);
                Vec3d e20 = pos(t[0]) - pos(t[2]);

                Vec3d gamma = -(e01 * (*m_Gamma)[t[2]].get(l) + e12 * (*m_Gamma)[t[0]].get(l)
                                + e20 * (*m_Gamma)[t[1]].get(l));

                Vec3d dx = x - xp;
                //                double dxn = dx.norm();
                double dxn = sqrt(dx.squaredNorm() + m_delta * m_delta);

                v += gamma.cross(dx) / (dxn * dxn * dxn);
                //                v += gamma.cross(dx) / (dxn * dxn * dxn) * (1 - exp(-dxn /
                //                m_delta));
            }

            v /= (4 * M_PI);
            m_dbg_v2[i] = v;
        }
    }

    // compute the mean curvatures
    if (true)
    {
        m_dbg_v1.clear();
        m_dbg_e1.clear();
        m_dbg_v1.resize(mesh().nv(), std::vector<double>(m_nregion, 0));
        m_dbg_e1.resize(mesh().ne(), std::vector<double>(m_nregion, 0));

        std::vector<std::vector<double>> curvature(
          mesh().ne(), std::vector<double>(m_nregion, 0)); // edge-aligned curvature (signed scalar)
        for (size_t region = 0; region < m_nregion; region++)
        {
            for (size_t i = 0; i < mesh().ne(); i++)
            {
                std::vector<size_t>
                  incident_faces; // faces incident to edge i that have the label of interest
                                  // (assume there are only two of them for now; this can be false
                                  // only when complex collision prevents immediate T1 resolution,
                                  // which is not expected to happen for bubble complexes.)
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
                //                assert(incident_faces.size() == 2);
                if (incident_faces.size() != 2)
                {
                    //                    std::cout << "Warning: incident_faces.size() != 2" <<
                    //                    std::endl; std::cout << "e = " << i << " region = " <<
                    //                    region << " n = " << incident_faces.size() << ": "; for
                    //                    (size_t k = 0; k < incident_faces.size(); k++) std::cout
                    //                    << mesh().m_edge_to_triangle_map[i][incident_faces[k]] <<
                    //                    " "; std::cout << std::endl;
                    curvature[i][region] = 0;
                    continue;
                }
                bool nonmanifold = (incident_regions.size() > 2);

                //                assert(mesh().m_edge_to_triangle_map[i].size() == 2);
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
                      std::swap(t0, t1); // the region of interest should be to the CCW direction of
                                         // t0 when looking along the direciton of edge i

                assert(mesh().get_triangle_label(ti0)[mesh().oriented(v0, v1, t0) ? 0 : 1]
                       == region);
                assert(mesh().get_triangle_label(ti1)[mesh().oriented(v0, v1, t1) ? 1 : 0]
                       == region);

                Vec3d n0 = (pos(t0[1]) - pos(t0[0])).cross(pos(t0[2]) - pos(t0[0]));
                if (mesh().get_triangle_label(ti0)[1] == region)
                    n0 = -n0; // n0 should point away from the region of interest

                Vec3d n1 = (pos(t1[1]) - pos(t1[0])).cross(pos(t1[2]) - pos(t1[0]));
                if (mesh().get_triangle_label(ti1)[1] == region)
                    n1 = -n1; // n1 should point away from the region of interest

                double edgeArea = (n0.norm() + n1.norm()) / 2 / 3;
                n0.normalize();
                n1.normalize();

                //                double curvature_i = 0;
                //                if (nonmanifold)
                //                {
                //                    curvature_i = n0.cross(n1).dot(et) / (et.norm() * m_delta); //
                //                    integral curvature along the curve orghotonal to the edge in
                //                    plane
                //                } else
                //                {
                //                    curvature_i = n0.cross(n1).dot(et) / edgeArea;
                ////                    curvature_i = angleAroundAxis(n0, n1, et) * et.norm() /
                /// edgeArea;
                //                }

                double curvature_i = n0.cross(n1).dot(et);
                //                double curvature_i = (n0 - n1).dot((n0 +
                //                n1).normalized().cross(et));

                curvature[i][region] = curvature_i;
                //                if (region == 1)
                //                    m_dbg_e1[i] = curvature_i;
                m_dbg_e1[i][region] = curvature_i;
            }

            //            for (size_t i = 0; i < mesh().nv(); i++)
            //            {
            //                double mean_curvature = 0;
            //                for (size_t j = 0; j < mesh().m_vertex_to_edge_map[i].size(); j++)
            //                    mean_curvature +=
            //                    curvature[mesh().m_vertex_to_edge_map[i][j]][region];
            //                mean_curvature /= mesh().m_vertex_to_edge_map[i].size();
            //
            //                (*m_Gamma)[i][region] += simOptions().sigma * mean_curvature * dt;
            //                m_dbg_v1[i][region] = mean_curvature;
            //            }

            for (size_t i = 0; i < mesh().nv(); i++)
            {
                double mean_curvature = 0;

                Mat3d second_fundamental_form = Mat3d::Zero();
                int counter = 0;
                double vertex_area = 0;
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
                        Vec3d et =
                          (pos(mesh().m_edges[e][1]) - pos(mesh().m_edges[e][0])).normalized();
                        second_fundamental_form += et * et.transpose() * curvature[e][region];
                        mean_curvature += curvature[e][region];
                        counter++;
                    }
                }

                for (size_t j = 0; j < mesh().m_vertex_to_triangle_map[i].size(); j++)
                {
                    const LosTopos::Vec3st& t =
                      mesh().m_tris[mesh().m_vertex_to_triangle_map[i][j]];
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

                //                if (counter == 0)
                //                    mean_curvature = 0;
                //                else
                ////                    mean_curvature = second_fundamental_form.trace() / counter /
                /// 3;
                //                    mean_curvature /= counter;

                mean_curvature = mean_curvature / (vertex_area * 2);

                m_dbg_v1[i][region] = mean_curvature;
            }
        }
    }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  Callbacks
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool
VS3D::generate_collapsed_position(LosTopos::SurfTrack& st,
                                  size_t v0,
                                  size_t v1,
                                  LosTopos::Vec3d& pos)
{
    if (st.vertex_is_any_solid(v0) && st.vertex_is_any_solid(v1))
    {
        return false;
    }
    else if (st.vertex_is_any_solid(v0))
    {
        pos = st.pm_positions[v0];
        return true;
    }
    else if (st.vertex_is_any_solid(v1))
    {
        pos = st.pm_positions[v1];
        return true;
    }
    else
    {
        pos = (st.pm_positions[v0] + st.pm_positions[v1]) / 2;
        return true;
    }
}

bool
VS3D::generate_split_position(LosTopos::SurfTrack& st, size_t v0, size_t v1, LosTopos::Vec3d& pos)
{
    std::cout << "solid callback: generate split position: " << v0 << " " << v1 << " "
              << (st.vertex_is_any_solid(v0) && st.vertex_is_any_solid(v1)) << std::endl;
    pos = (st.pm_positions[v0] + st.pm_positions[v1]) / 2;
    if (st.vertex_is_any_solid(v0) && st.vertex_is_any_solid(v1))
        return false;
    else
        return true;
}

LosTopos::Vec3c
VS3D::generate_collapsed_solid_label(LosTopos::SurfTrack& st,
                                     size_t v0,
                                     size_t v1,
                                     const LosTopos::Vec3c& label0,
                                     const LosTopos::Vec3c& label1)
{
    return LosTopos::Vec3c(1, 1, 1);
}

LosTopos::Vec3c
VS3D::generate_split_solid_label(LosTopos::SurfTrack& st,
                                 size_t v0,
                                 size_t v1,
                                 const LosTopos::Vec3c& label0,
                                 const LosTopos::Vec3c& label1)
{
    if (st.vertex_is_any_solid(v0) && st.vertex_is_any_solid(v1))
        return LosTopos::Vec3c(1, 1, 1);
    else if (st.vertex_is_any_solid(v0))
        return LosTopos::Vec3c(0, 0, 0);
    else if (st.vertex_is_any_solid(v1))
        return LosTopos::Vec3c(0, 0, 0);
    else
        return LosTopos::Vec3c(0, 0, 0);
}

bool
VS3D::generate_edge_popped_positions(LosTopos::SurfTrack& st,
                                     size_t oldv,
                                     const LosTopos::Vec2i& cut,
                                     LosTopos::Vec3d& pos_upper,
                                     LosTopos::Vec3d& pos_lower)
{
    return false;
}

bool
VS3D::generate_vertex_popped_positions(LosTopos::SurfTrack& st,
                                       size_t oldv,
                                       int A,
                                       int B,
                                       LosTopos::Vec3d& pos_a,
                                       LosTopos::Vec3d& pos_b)
{
    return false;
}

bool
VS3D::solid_edge_is_feature(const LosTopos::SurfTrack& st, size_t e)
{
    return false;
}

LosTopos::Vec3d
VS3D::sampleVelocity(LosTopos::Vec3d& pos)
{
    return LosTopos::Vec3d(0, 0, 0);
}

bool
VS3D::sampleDirectionalDivergence(const LosTopos::Vec3d& pos,
                                  const LosTopos::Vec3d& dir,
                                  double& output)
{
    return false;
}

struct CollapseTempData
{
    size_t v0;
    size_t v1;

    Vec3d old_x0;
    Vec3d old_x1;

    Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> v0_incident_region_pairs;
    Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> v1_incident_region_pairs;
};

void
VS3D::pre_collapse(const LosTopos::SurfTrack& st, size_t e, void** data)
{
    CollapseTempData* td = new CollapseTempData;
    td->v0 = st.m_mesh.m_edges[e][0];
    td->v1 = st.m_mesh.m_edges[e][1];

    td->old_x0 = vc(st.pm_positions[td->v0]);
    td->old_x1 = vc(st.pm_positions[td->v1]);

    td->v0_incident_region_pairs.setZero(m_nregion, m_nregion);
    for (size_t i = 0; i < st.m_mesh.m_vertex_to_triangle_map[td->v0].size(); i++)
    {
        LosTopos::Vec2i l =
          st.m_mesh.get_triangle_label(st.m_mesh.m_vertex_to_triangle_map[td->v0][i]);
        td->v0_incident_region_pairs(l[0], l[1]) = td->v0_incident_region_pairs(l[1], l[0]) = true;
    }
    td->v1_incident_region_pairs.setZero(m_nregion, m_nregion);
    for (size_t i = 0; i < st.m_mesh.m_vertex_to_triangle_map[td->v1].size(); i++)
    {
        LosTopos::Vec2i l =
          st.m_mesh.get_triangle_label(st.m_mesh.m_vertex_to_triangle_map[td->v1][i]);
        td->v1_incident_region_pairs(l[0], l[1]) = td->v1_incident_region_pairs(l[1], l[0]) = true;
    }

    *data = (void*)td;
    std::cout << "pre collapse: " << e << ": " << td->v0 << " " << td->v1 << std::endl;
}

void
VS3D::post_collapse(const LosTopos::SurfTrack& st, size_t e, size_t merged_vertex, void* data)
{
    CollapseTempData* td = (CollapseTempData*)data;
    std::cout << "post collapse: " << e << ": " << td->v0 << " " << td->v1 << " => "
              << merged_vertex << std::endl;
    assert((st.m_mesh.vertex_is_deleted(td->v0) && merged_vertex == td->v1)
           || (st.m_mesh.vertex_is_deleted(td->v1) && merged_vertex == td->v0));

    Vec3d merged_x = vc(st.pm_positions[merged_vertex]);
    double s = (merged_x - td->old_x0).dot(td->old_x1 - td->old_x0)
               / (td->old_x1 - td->old_x0).squaredNorm();

    Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> merged_vertex_incident_region_pairs;
    merged_vertex_incident_region_pairs.setZero(m_nregion, m_nregion);
    for (size_t i = 0; i < st.m_mesh.m_vertex_to_triangle_map[merged_vertex].size(); i++)
    {
        LosTopos::Vec2i l =
          st.m_mesh.get_triangle_label(st.m_mesh.m_vertex_to_triangle_map[merged_vertex][i]);
        merged_vertex_incident_region_pairs(l[0], l[1]) =
          merged_vertex_incident_region_pairs(l[1], l[0]) = true;
    }

    // for region pairs not existing in original v0 and v1, we cannot simply assume 0 circulation
    // because the neighbors that do have those region pairs
    //  may have accumulated some amount of circulations, and filling in 0 will cause vorticity
    //  spikes. The information that should be filled in here can only come from neighbors.
    GammaType newGamma(m_nregion);
    newGamma.setZero();
    for (int i = 0; i < m_nregion; i++)
    {
        for (int j = i + 1; j < m_nregion; j++)
        {
            if (merged_vertex_incident_region_pairs(i, j))
            {
                if (!td->v0_incident_region_pairs(i, j) && !td->v1_incident_region_pairs(i, j))
                {
                    double neighborhood_mean = 0;
                    int neighborhood_counter = 0;
                    for (size_t k = 0; k < st.m_mesh.m_vertex_to_edge_map[merged_vertex].size();
                         k++)
                    {
                        LosTopos::Vec2st e =
                          st.m_mesh.m_edges[st.m_mesh.m_vertex_to_edge_map[merged_vertex][k]];
                        size_t vother = (e[0] == merged_vertex ? e[1] : e[0]);
                        bool incident_to_this_region_pair = false;
                        for (size_t l = 0; l < st.m_mesh.m_vertex_to_triangle_map[vother].size();
                             l++)
                        {
                            LosTopos::Vec2i ll = st.m_mesh.get_triangle_label(
                              st.m_mesh.m_vertex_to_triangle_map[vother][l]);
                            if ((ll[0] == i && ll[1] == j) || (ll[0] == j && ll[1] == i))
                            {
                                incident_to_this_region_pair = true;
                                break;
                            }
                        }

                        if (incident_to_this_region_pair)
                        {
                            neighborhood_mean += (*m_Gamma)[vother].get(i, j);
                            neighborhood_counter++;
                        }
                    }
                    if (neighborhood_counter != 0)
                        neighborhood_mean /= neighborhood_counter;

                    newGamma.set(i, j, neighborhood_mean);
                }
                else if (!td->v0_incident_region_pairs(i, j))
                {
                    newGamma.set(i, j, (*m_Gamma)[td->v1].get(i, j));
                }
                else if (!td->v1_incident_region_pairs(i, j))
                {
                    newGamma.set(i, j, (*m_Gamma)[td->v0].get(i, j));
                }
                else
                {
                    newGamma.set(i,
                                 j,
                                 (*m_Gamma)[td->v0].get(i, j) * (1 - s)
                                   + (*m_Gamma)[td->v1].get(i, j) * s);
                }
            }
        }
    }
    (*m_Gamma)[merged_vertex] = newGamma;
}

struct SplitTempData
{
    size_t v0;
    size_t v1;

    Vec3d old_x0;
    Vec3d old_x1;

    Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> v0_incident_region_pairs;
    Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> v1_incident_region_pairs;
};

void
VS3D::pre_split(const LosTopos::SurfTrack& st, size_t e, void** data)
{
    SplitTempData* td = new SplitTempData;
    td->v0 = st.m_mesh.m_edges[e][0];
    td->v1 = st.m_mesh.m_edges[e][1];

    td->old_x0 = vc(st.pm_positions[td->v0]);
    td->old_x1 = vc(st.pm_positions[td->v1]);

    td->v0_incident_region_pairs.setZero(m_nregion, m_nregion);
    for (size_t i = 0; i < st.m_mesh.m_vertex_to_triangle_map[td->v0].size(); i++)
    {
        LosTopos::Vec2i l =
          st.m_mesh.get_triangle_label(st.m_mesh.m_vertex_to_triangle_map[td->v0][i]);
        td->v0_incident_region_pairs(l[0], l[1]) = td->v0_incident_region_pairs(l[1], l[0]) = true;
    }
    td->v1_incident_region_pairs.setZero(m_nregion, m_nregion);
    for (size_t i = 0; i < st.m_mesh.m_vertex_to_triangle_map[td->v1].size(); i++)
    {
        LosTopos::Vec2i l =
          st.m_mesh.get_triangle_label(st.m_mesh.m_vertex_to_triangle_map[td->v1][i]);
        td->v1_incident_region_pairs(l[0], l[1]) = td->v1_incident_region_pairs(l[1], l[0]) = true;
    }

    *data = (void*)td;
    std::cout << "pre split: " << e << ": " << td->v0 << " " << td->v1 << std::endl;
}

void
VS3D::post_split(const LosTopos::SurfTrack& st, size_t e, size_t new_vertex, void* data)
{
    SplitTempData* td = (SplitTempData*)data;
    std::cout << "post split: " << e << ": " << td->v0 << " " << td->v1 << " => " << new_vertex
              << std::endl;

    Vec3d midpoint_x = vc(st.pm_positions[new_vertex]);
    double s = (midpoint_x - td->old_x0).dot(td->old_x1 - td->old_x0)
               / (td->old_x1 - td->old_x0).squaredNorm();

    Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> new_vertex_incident_region_pairs;
    new_vertex_incident_region_pairs.setZero(m_nregion, m_nregion);
    for (size_t i = 0; i < st.m_mesh.m_vertex_to_triangle_map[new_vertex].size(); i++)
    {
        LosTopos::Vec2i l =
          st.m_mesh.get_triangle_label(st.m_mesh.m_vertex_to_triangle_map[new_vertex][i]);
        new_vertex_incident_region_pairs(l[0], l[1]) =
          new_vertex_incident_region_pairs(l[1], l[0]) = true;
    }

    // for region pairs not existing in original v0 and v1, we cannot simply assume 0 circulation
    // because the neighbors that do have those region pairs
    //  may have accumulated some amount of circulations, and filling in 0 will cause vorticity
    //  spikes. The information that should be filled in here can only come from neighbors.
    GammaType newGamma(m_nregion);
    newGamma.setZero();
    for (int i = 0; i < m_nregion; i++)
    {
        for (int j = i + 1; j < m_nregion; j++)
        {
            if (new_vertex_incident_region_pairs(i, j))
            {
                if (!td->v0_incident_region_pairs(i, j) && !td->v1_incident_region_pairs(i, j))
                {
                    double neighborhood_mean = 0;
                    int neighborhood_counter = 0;
                    for (size_t k = 0; k < st.m_mesh.m_vertex_to_edge_map[new_vertex].size(); k++)
                    {
                        LosTopos::Vec2st e =
                          st.m_mesh.m_edges[st.m_mesh.m_vertex_to_edge_map[new_vertex][k]];
                        size_t vother = (e[0] == new_vertex ? e[1] : e[0]);

                        bool incident_to_this_region_pair = false;
                        for (size_t l = 0; l < st.m_mesh.m_vertex_to_triangle_map[vother].size();
                             l++)
                        {
                            LosTopos::Vec2i ll = st.m_mesh.get_triangle_label(
                              st.m_mesh.m_vertex_to_triangle_map[vother][l]);
                            if ((ll[0] == i && ll[1] == j) || (ll[0] == j && ll[1] == i))
                            {
                                incident_to_this_region_pair = true;
                                break;
                            }
                        }

                        if (incident_to_this_region_pair)
                        {
                            neighborhood_mean += (*m_Gamma)[vother].get(i, j);
                            neighborhood_counter++;
                        }
                    }
                    if (neighborhood_counter != 0)
                        neighborhood_mean /= neighborhood_counter;

                    newGamma.set(i, j, neighborhood_mean);
                }
                else if (!td->v0_incident_region_pairs(i, j))
                {
                    newGamma.set(i, j, (*m_Gamma)[td->v1].get(i, j));
                }
                else if (!td->v1_incident_region_pairs(i, j))
                {
                    newGamma.set(i, j, (*m_Gamma)[td->v0].get(i, j));
                }
                else
                {
                    newGamma.set(i,
                                 j,
                                 (*m_Gamma)[td->v0].get(i, j) * (1 - s)
                                   + (*m_Gamma)[td->v1].get(i, j) * s);
                }
            }
        }
    }
    (*m_Gamma)[new_vertex] = newGamma;
}

void
VS3D::pre_flip(const LosTopos::SurfTrack& st, size_t e, void** data)
{
}

void
VS3D::post_flip(const LosTopos::SurfTrack& st, size_t e, void* data)
{
}

struct T1TempData
{
    std::vector<size_t> neighbor_verts;
    std::vector<Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic>> neighbor_region_pairs;
};

void
VS3D::pre_t1(const LosTopos::SurfTrack& st, size_t v, void** data)
{
    T1TempData* td = new T1TempData;

    for (size_t i = 0; i < st.m_mesh.m_vertex_to_edge_map[v].size(); i++)
    {
        LosTopos::Vec2st e = st.m_mesh.m_edges[st.m_mesh.m_vertex_to_edge_map[v][i]];
        size_t vother = (e[0] == v ? e[1] : e[0]);
        td->neighbor_verts.push_back(vother);

        Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> rps;
        rps.setZero(m_nregion, m_nregion);
        for (size_t j = 0; j < st.m_mesh.m_vertex_to_triangle_map[vother].size(); j++)
        {
            LosTopos::Vec2i l =
              st.m_mesh.get_triangle_label(st.m_mesh.m_vertex_to_triangle_map[vother][j]);
            rps(l[0], l[1]) = rps(l[1], l[0]) = true;
        }
        td->neighbor_region_pairs.push_back(rps);
    }

    *data = (void*)td;
}

void
VS3D::post_t1(const LosTopos::SurfTrack& st, size_t v, size_t a, size_t b, void* data)
{
    std::cout << "v = " << v << " -> " << a << " " << b << std::endl;
    T1TempData* td = (T1TempData*)data;

    (*m_Gamma)[a] = (*m_Gamma)[v];
    (*m_Gamma)[b] = (*m_Gamma)[v];

    std::cout << "v Gammas: " << std::endl << (*m_Gamma)[v].values << std::endl;

    for (size_t i = 0; i < td->neighbor_verts.size(); i++)
    {
        size_t vother = td->neighbor_verts[i];

        Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> rps;
        rps.setZero(m_nregion, m_nregion);
        for (size_t j = 0; j < st.m_mesh.m_vertex_to_triangle_map[vother].size(); j++)
        {
            LosTopos::Vec2i l =
              st.m_mesh.get_triangle_label(st.m_mesh.m_vertex_to_triangle_map[vother][j]);
            rps(l[0], l[1]) = rps(l[1], l[0]) = true;
        }

        std::cout << vother << std::endl
                  << td->neighbor_region_pairs[i] << std::endl
                  << "-> " << std::endl
                  << rps << std::endl;

        for (int j = 0; j < m_nregion; j++)
            for (int k = j + 1; k < m_nregion; k++)
                if (rps(j, k) && !td->neighbor_region_pairs[i](j, k))
                    (*m_Gamma)[vother].set(j, k, (*m_Gamma)[v].get(j, k));

        std::cout << "Gammas: " << std::endl << (*m_Gamma)[vother].values << std::endl;
    }
}

struct FaceSplitTempData
{
    size_t v0;
    size_t v1;
    size_t v2;
};

void
VS3D::pre_facesplit(const LosTopos::SurfTrack& st, size_t f, void** data)
{
    FaceSplitTempData* td = new FaceSplitTempData;

    td->v0 = st.m_mesh.get_triangle(f)[0];
    td->v1 = st.m_mesh.get_triangle(f)[1];
    td->v2 = st.m_mesh.get_triangle(f)[2];

    *data = (void*)td;
}

void
VS3D::post_facesplit(const LosTopos::SurfTrack& st, size_t f, size_t new_vertex, void* data)
{
    FaceSplitTempData* td = (FaceSplitTempData*)data;

    (*m_Gamma)[new_vertex] = GammaType(m_nregion);
    (*m_Gamma)[new_vertex].values =
      ((*m_Gamma)[td->v0].values + (*m_Gamma)[td->v1].values + (*m_Gamma)[td->v2].values) / 3;
}

struct SnapTempData
{
    size_t v0;
    size_t v1;

    Vec3d old_x0;
    Vec3d old_x1;

    Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> v0_incident_region_pairs;
    Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> v1_incident_region_pairs;
};

void
VS3D::pre_snap(const LosTopos::SurfTrack& st, size_t v0, size_t v1, void** data)
{
    SnapTempData* td = new SnapTempData;
    td->v0 = v0;
    td->v1 = v1;

    td->old_x0 = vc(st.pm_positions[v0]);
    td->old_x1 = vc(st.pm_positions[v1]);

    td->v0_incident_region_pairs.setZero(m_nregion, m_nregion);
    for (size_t i = 0; i < st.m_mesh.m_vertex_to_triangle_map[td->v0].size(); i++)
    {
        LosTopos::Vec2i l =
          st.m_mesh.get_triangle_label(st.m_mesh.m_vertex_to_triangle_map[td->v0][i]);
        td->v0_incident_region_pairs(l[0], l[1]) = td->v0_incident_region_pairs(l[1], l[0]) = true;
    }
    td->v1_incident_region_pairs.setZero(m_nregion, m_nregion);
    for (size_t i = 0; i < st.m_mesh.m_vertex_to_triangle_map[td->v1].size(); i++)
    {
        LosTopos::Vec2i l =
          st.m_mesh.get_triangle_label(st.m_mesh.m_vertex_to_triangle_map[td->v1][i]);
        td->v1_incident_region_pairs(l[0], l[1]) = td->v1_incident_region_pairs(l[1], l[0]) = true;
    }

    *data = (void*)td;
    std::cout << "pre snap: " << v0 << " " << v1 << std::endl;
}

void
VS3D::post_snap(const LosTopos::SurfTrack& st, size_t v_kept, size_t v_deleted, void* data)
{
    SnapTempData* td = (SnapTempData*)data;
    std::cout << "post snap: " << td->v0 << " " << td->v1 << " => " << v_kept << std::endl;
    assert((td->v0 == v_kept && td->v1 == v_deleted) || (td->v1 == v_kept && td->v0 == v_deleted));
    assert(v_kept != v_deleted);
    assert(st.m_mesh.vertex_is_deleted(v_deleted));
    assert(!st.m_mesh.vertex_is_deleted(v_kept));

    Vec3d merged_x = vc(st.pm_positions[v_kept]);
    double s = (merged_x - td->old_x0).dot(td->old_x1 - td->old_x0)
               / (td->old_x1 - td->old_x0).squaredNorm();

    Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> merged_vertex_incident_region_pairs;
    merged_vertex_incident_region_pairs.setZero(m_nregion, m_nregion);
    for (size_t i = 0; i < st.m_mesh.m_vertex_to_triangle_map[v_kept].size(); i++)
    {
        LosTopos::Vec2i l =
          st.m_mesh.get_triangle_label(st.m_mesh.m_vertex_to_triangle_map[v_kept][i]);
        std::cout << "triangle " << st.m_mesh.m_vertex_to_triangle_map[v_kept][i]
                  << " label = " << l << std::endl;
        merged_vertex_incident_region_pairs(l[0], l[1]) =
          merged_vertex_incident_region_pairs(l[1], l[0]) = true;
    }

    std::cout << "v0 incident region pairs: " << std::endl
              << td->v0_incident_region_pairs << std::endl;
    std::cout << "v1 incident region pairs: " << std::endl
              << td->v1_incident_region_pairs << std::endl;
    std::cout << "merged vertex incident region pairs: " << std::endl
              << merged_vertex_incident_region_pairs << std::endl;

    // for region pairs not existing in original v0 and v1, we cannot simply assume 0 circulation
    // because the neighbors that do have those region pairs
    //  may have accumulated some amount of circulations, and filling in 0 will cause vorticity
    //  spikes. The information that should be filled in here can only come from neighbors.
    GammaType newGamma(m_nregion);
    newGamma.setZero();
    for (int i = 0; i < m_nregion; i++)
    {
        for (int j = i + 1; j < m_nregion; j++)
        {
            if (merged_vertex_incident_region_pairs(i, j))
            {
                if (!td->v0_incident_region_pairs(i, j) && !td->v1_incident_region_pairs(i, j))
                {
                    std::cout << "region pair " << i << " " << j
                              << " is being computed from 1-ring neighbors" << std::endl;
                    double neighborhood_mean = 0;
                    int neighborhood_counter = 0;
                    for (size_t k = 0; k < st.m_mesh.m_vertex_to_edge_map[v_kept].size(); k++)
                    {
                        LosTopos::Vec2st e =
                          st.m_mesh.m_edges[st.m_mesh.m_vertex_to_edge_map[v_kept][k]];
                        size_t vother = (e[0] == v_kept ? e[1] : e[0]);

                        bool incident_to_this_region_pair = false;
                        for (size_t l = 0; l < st.m_mesh.m_vertex_to_triangle_map[vother].size();
                             l++)
                        {
                            LosTopos::Vec2i ll = st.m_mesh.get_triangle_label(
                              st.m_mesh.m_vertex_to_triangle_map[vother][l]);
                            if ((ll[0] == i && ll[1] == j) || (ll[0] == j && ll[1] == i))
                            {
                                incident_to_this_region_pair = true;
                                break;
                            }
                        }

                        std::cout << "vother = " << vother
                                  << " incident = " << incident_to_this_region_pair << std::endl;

                        if (incident_to_this_region_pair)
                        {
                            neighborhood_mean += (*m_Gamma)[vother].get(i, j);
                            neighborhood_counter++;
                        }
                    }
                    if (neighborhood_counter != 0)
                        neighborhood_mean /= neighborhood_counter;

                    newGamma.set(i, j, neighborhood_mean);
                }
                else if (!td->v0_incident_region_pairs(i, j))
                {
                    newGamma.set(i, j, (*m_Gamma)[td->v1].get(i, j));
                }
                else if (!td->v1_incident_region_pairs(i, j))
                {
                    newGamma.set(i, j, (*m_Gamma)[td->v0].get(i, j));
                }
                else
                {
                    newGamma.set(i,
                                 j,
                                 (*m_Gamma)[td->v0].get(i, j) * (1 - s)
                                   + (*m_Gamma)[td->v1].get(i, j) * s);
                }
            }
        }
    }
    std::cout << "v0 Gamma = " << std::endl << (*m_Gamma)[td->v0].values << std::endl;
    std::cout << "v1 Gamma = " << std::endl << (*m_Gamma)[td->v1].values << std::endl;
    std::cout << "new Gamma = " << std::endl << newGamma.values << std::endl;
    (*m_Gamma)[v_kept] = newGamma;
}

void
VS3D::accumulateGradU(VectorXs& F, const VectorXs& dx, const VectorXs& dv)
{
    const int ndof = m_constrained_positions.size() * 3;

    // Accumulate all energy gradients
    if (dx.size() == 0)
    {
        for (std::vector<Force*>::size_type i = 0; i < m_forces.size(); ++i)
            m_forces[i]->addGradEToTotal(
              Eigen::Map<VectorXs>((double*)&m_constrained_positions[0], ndof),
              Eigen::Map<VectorXs>((double*)&m_constrained_velocities[0], ndof),
              Eigen::Map<VectorXs>((double*)&m_constrained_mass[0], ndof),
              F);
    }
    else
    {
        VectorXs nx = Eigen::Map<VectorXs>((double*)&m_constrained_positions[0], ndof) + dx;
        VectorXs nv = Eigen::Map<VectorXs>((double*)&m_constrained_velocities[0], ndof) + dv;
        for (std::vector<Force*>::size_type i = 0; i < m_forces.size(); ++i)
            m_forces[i]->addGradEToTotal(
              nx, nv, Eigen::Map<VectorXs>((double*)&m_constrained_mass[0], ndof), F);
    }
}

void
VS3D::accumulateddUdxdx(TripletXs& A, const VectorXs& dx, const VectorXs& dv)
{
    const int ndof = m_constrained_positions.size() * 3;
    if (dx.size() == 0)
    {
        for (std::vector<Force*>::size_type i = 0; i < m_forces.size(); ++i)
            m_forces[i]->addHessXToTotal(
              Eigen::Map<VectorXs>((double*)&m_constrained_positions[0], ndof),
              Eigen::Map<VectorXs>((double*)&m_constrained_velocities[0], ndof),
              Eigen::Map<VectorXs>((double*)&m_constrained_mass[0], ndof),
              A);
    }
    else
    {
        VectorXs nx = Eigen::Map<VectorXs>((double*)&m_constrained_positions[0], ndof) + dx;
        VectorXs nv = Eigen::Map<VectorXs>((double*)&m_constrained_velocities[0], ndof) + dv;
        for (std::vector<Force*>::size_type i = 0; i < m_forces.size(); ++i)
            m_forces[i]->addHessXToTotal(
              nx, nv, Eigen::Map<VectorXs>((double*)&m_constrained_mass[0], ndof), A);
    }
}

// Kind of a misnomer.
void
VS3D::accumulateddUdxdv(TripletXs& A, const VectorXs& dx, const VectorXs& dv)
{
    const int ndof = m_constrained_positions.size() * 3;
    if (dx.size() == 0)
    {
        for (std::vector<Force*>::size_type i = 0; i < m_forces.size(); ++i)
            m_forces[i]->addHessVToTotal(
              Eigen::Map<VectorXs>((double*)&m_constrained_positions[0], ndof),
              Eigen::Map<VectorXs>((double*)&m_constrained_velocities[0], ndof),
              Eigen::Map<VectorXs>((double*)&m_constrained_mass[0], ndof),
              A);
    }
    else
    {
        VectorXs nx = Eigen::Map<VectorXs>((double*)&m_constrained_positions[0], ndof) + dx;
        VectorXs nv = Eigen::Map<VectorXs>((double*)&m_constrained_velocities[0], ndof) + dv;
        for (std::vector<Force*>::size_type i = 0; i < m_forces.size(); ++i)
            m_forces[i]->addHessVToTotal(
              nx, nv, Eigen::Map<VectorXs>((double*)&m_constrained_mass[0], ndof), A);
    }
}

void
VS3D::preCompute(const VectorXs& dx, const VectorXs& dv, const scalar& dt)
{
    const int ndof = m_constrained_positions.size() * 3;
    if (dx.size() == 0)
    {
        for (std::vector<Force*>::size_type i = 0; i < m_forces.size(); ++i)
            m_forces[i]->preCompute(
              Eigen::Map<VectorXs>((double*)&m_constrained_positions[0], ndof),
              Eigen::Map<VectorXs>((double*)&m_constrained_velocities[0], ndof),
              Eigen::Map<VectorXs>((double*)&m_constrained_mass[0], ndof),
              dt);
    }
    else
    {
        VectorXs nx = Eigen::Map<VectorXs>((double*)&m_constrained_positions[0], ndof) + dx;
        VectorXs nv = Eigen::Map<VectorXs>((double*)&m_constrained_velocities[0], ndof) + dv;
        for (std::vector<Force*>::size_type i = 0; i < m_forces.size(); ++i)
            m_forces[i]->preCompute(
              nx, nv, Eigen::Map<VectorXs>((double*)&m_constrained_mass[0], ndof), dt);
    }
}

void
VS3D::pre_smoothing(const LosTopos::SurfTrack& st, void** data)
{
}

void
VS3D::post_smoothing(const LosTopos::SurfTrack& st, void* data)
{
}
