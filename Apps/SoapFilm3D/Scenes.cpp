//
//  Scenes.cpp
//  MultiTracker
//
//  Created by Fang Da on 15/1/26.
//
//

#include "Scenes.h"
#include "MeshIO.h"
#include "Sim.h"
#include "SimOptions.h"
#include "VS3D.h"
#include <boost/multi_array.hpp>
#include <map>
#include <random>

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  Scene-specific initialization
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
namespace
{
class OrderComp
{
  public:
    bool operator()(const std::pair<int, double>& o1, const std::pair<int, double>& o2) const
    {
        return o1.second < o2.second;
    }
};

// subdivide every triangle in the mesh into four triangles, by introducing three new vertices at
// edge midpoints after subdivision, vs will be expanded (retaining all original vertices), while fs
// and ls will be replaced (no original faces remain) if r is positive, new vertices will be
// projected onto the sphere centered at center, with a radius of r (for ICO sphere generation)
void
subdivide(const Vec3d& center,
          double r,
          std::vector<Vec3d>& vs,
          std::vector<Vec3i>& fs,
          std::vector<Vec2i>& ls)
{
    std::vector<Vec3i> new_fs;
    std::vector<Vec2i> new_ls;
    std::map<std::pair<int, int>, int> new_vs_map;

    // map edge-midpoint to new vertices
    for (size_t j = 0; j < fs.size(); j++)
    {
        Vec3i& f = fs[j];
        int v0 = f[0];
        int v1 = f[1];
        int v2 = f[2];

        std::pair<int, int> p01 =
          (v0 < v1 ? std::pair<int, int>(v0, v1) : std::pair<int, int>(v1, v0));
        std::pair<int, int> p12 =
          (v1 < v2 ? std::pair<int, int>(v1, v2) : std::pair<int, int>(v2, v1));
        std::pair<int, int> p20 =
          (v2 < v0 ? std::pair<int, int>(v2, v0) : std::pair<int, int>(v0, v2));

        new_vs_map[p01] = 0;
        new_vs_map[p12] = 0;
        new_vs_map[p20] = 0;
    }

    // create the new vertices
    for (std::map<std::pair<int, int>, int>::iterator j = new_vs_map.begin(); j != new_vs_map.end();
         j++)
    {
        j->second = vs.size();
        if (r > 0)
            vs.push_back(((vs[j->first.first] + vs[j->first.second]) / 2 - center).normalized() * r
                         + center);
        else
            vs.push_back((vs[j->first.first] + vs[j->first.second]) / 2);
    }

    // triangulate
    for (size_t j = 0; j < fs.size(); j++)
    {
        Vec3i& f = fs[j];
        int v0 = f[0];
        int v1 = f[1];
        int v2 = f[2];

        std::pair<int, int> p01 =
          (v0 < v1 ? std::pair<int, int>(v0, v1) : std::pair<int, int>(v1, v0));
        std::pair<int, int> p12 =
          (v1 < v2 ? std::pair<int, int>(v1, v2) : std::pair<int, int>(v2, v1));
        std::pair<int, int> p20 =
          (v2 < v0 ? std::pair<int, int>(v2, v0) : std::pair<int, int>(v0, v2));
        int nv01 = new_vs_map[p01];
        int nv12 = new_vs_map[p12];
        int nv20 = new_vs_map[p20];

        Vec2i& l = ls[j];
        new_fs.push_back(Vec3i(v0, nv01, nv20));
        new_ls.push_back(l);
        new_fs.push_back(Vec3i(nv01, v1, nv12));
        new_ls.push_back(l);
        new_fs.push_back(Vec3i(nv20, nv12, v2));
        new_ls.push_back(l);
        new_fs.push_back(Vec3i(nv12, nv20, nv01));
        new_ls.push_back(l);
    }

    fs = new_fs;
    ls = new_ls;
}

void
createIcoSphere(const Vec3d& center,
                double r,
                int subdivision,
                std::vector<Vec3d>& vs,
                std::vector<Vec3i>& fs,
                std::vector<Vec2i>& ls,
                const Vec2i& label = Vec2i(1, 0))
{
    // create the initial icosahedron
    double phi = (1.0 + sqrt(5.0)) / 2.0;
    double len = Vec3d(0, 1, phi).norm();

    vs.resize(12);
    vs[0] = center + r * Vec3d(0, 1, phi) / len;
    vs[1] = center + r * Vec3d(0, -1, phi) / len;
    vs[2] = center + r * Vec3d(0, 1, -phi) / len;
    vs[3] = center + r * Vec3d(0, -1, -phi) / len;
    vs[4] = center + r * Vec3d(1, phi, 0) / len;
    vs[5] = center + r * Vec3d(-1, phi, 0) / len;
    vs[6] = center + r * Vec3d(1, -phi, 0) / len;
    vs[7] = center + r * Vec3d(-1, -phi, 0) / len;
    vs[8] = center + r * Vec3d(phi, 0, 1) / len;
    vs[9] = center + r * Vec3d(phi, 0, -1) / len;
    vs[10] = center + r * Vec3d(-phi, 0, 1) / len;
    vs[11] = center + r * Vec3d(-phi, 0, -1) / len;

    fs.push_back(Vec3i(0, 1, 8));
    ls.push_back(label);
    fs.push_back(Vec3i(1, 0, 10));
    ls.push_back(label);
    fs.push_back(Vec3i(2, 3, 11));
    ls.push_back(label);
    fs.push_back(Vec3i(3, 2, 9));
    ls.push_back(label);

    fs.push_back(Vec3i(4, 5, 0));
    ls.push_back(label);
    fs.push_back(Vec3i(5, 4, 2));
    ls.push_back(label);
    fs.push_back(Vec3i(6, 7, 3));
    ls.push_back(label);
    fs.push_back(Vec3i(7, 6, 1));
    ls.push_back(label);

    fs.push_back(Vec3i(8, 9, 4));
    ls.push_back(label);
    fs.push_back(Vec3i(9, 8, 6));
    ls.push_back(label);
    fs.push_back(Vec3i(10, 11, 7));
    ls.push_back(label);
    fs.push_back(Vec3i(11, 10, 5));
    ls.push_back(label);

    fs.push_back(Vec3i(0, 8, 4));
    ls.push_back(label);
    fs.push_back(Vec3i(1, 6, 8));
    ls.push_back(label);
    fs.push_back(Vec3i(0, 5, 10));
    ls.push_back(label);
    fs.push_back(Vec3i(1, 10, 7));
    ls.push_back(label);

    fs.push_back(Vec3i(11, 3, 7));
    ls.push_back(label);
    fs.push_back(Vec3i(5, 2, 11));
    ls.push_back(label);
    fs.push_back(Vec3i(6, 3, 9));
    ls.push_back(label);
    fs.push_back(Vec3i(9, 2, 4));
    ls.push_back(label);

    for (int i = 0; i < subdivision; i++)
        subdivide(center, r, vs, fs, ls);
}

void
createCube(const Vec3d& position,
           double size,
           std::size_t number_subdivision,
           std::vector<Vec3d>& vs,
           std::vector<Vec3i>& fs,
           std::vector<Vec2i>& ls,
           const Vec2i& label)
{
    vs.clear();
    vs.reserve(8u);
    vs.push_back(position + size * Vec3d(0, 0, 0));
    vs.push_back(position + size * Vec3d(1, 0, 0));
    vs.push_back(position + size * Vec3d(1, 1, 0));
    vs.push_back(position + size * Vec3d(0, 1, 0));
    vs.push_back(position + size * Vec3d(0, 0, 1));
    vs.push_back(position + size * Vec3d(1, 0, 1));
    vs.push_back(position + size * Vec3d(1, 1, 1));
    vs.push_back(position + size * Vec3d(0, 1, 1));

    fs.clear();
    ls.clear();
    fs.reserve(12u);
    ls.reserve(12u);
    fs.push_back(Vec3i(0, 3, 2));
    ls.push_back(label);
    fs.push_back(Vec3i(2, 1, 0));
    ls.push_back(label);
    fs.push_back(Vec3i(0, 1, 5));
    ls.push_back(label);
    fs.push_back(Vec3i(4, 0, 5));
    ls.push_back(label);
    fs.push_back(Vec3i(1, 2, 6));
    ls.push_back(label);
    fs.push_back(Vec3i(6, 5, 1));
    ls.push_back(label);
    fs.push_back(Vec3i(2, 3, 7));
    ls.push_back(label);
    fs.push_back(Vec3i(7, 6, 2));
    ls.push_back(label);
    fs.push_back(Vec3i(0, 4, 3));
    ls.push_back(label);
    fs.push_back(Vec3i(4, 7, 3));
    ls.push_back(label);
    fs.push_back(Vec3i(4, 5, 6));
    ls.push_back(label);
    fs.push_back(Vec3i(7, 4, 6));
    ls.push_back(label);

    for (std::size_t i = 0; i < number_subdivision; ++i)
    {
        subdivide(Vec3d(0, 0, 0), 0, vs, fs, ls);
    }
}

void
addMesh(const std::vector<Vec3d>& mesh_vertices,
        const std::vector<Vec3i>& mesh_faces,
        const std::vector<Vec2i>& mesh_labels,
        std::vector<Vec3d>& vertices,
        std::vector<Vec3i>& faces,
        std::vector<Vec2i>& labels)
{
    Vec3i offset = vertices.size() * Vec3i::Ones();

    vertices.insert(vertices.end(), mesh_vertices.begin(), mesh_vertices.end());
    labels.insert(labels.end(), mesh_labels.begin(), mesh_labels.end());
    for (const Vec3i& mesh_face : mesh_faces)
    {
        faces.push_back(offset + mesh_face);
    }
}

template<typename SphereRange>
double
getSphereDistanceToSpheres(const Vec3d& center, double radius, const SphereRange& other_spheres)
{
    double min_distance = std::numeric_limits<double>::infinity();

    double distance;
    for (const auto& [other_center, other_radius] : other_spheres)
    {
        distance = (center - other_center).norm() - (other_radius + radius);
        min_distance = std::min(distance, min_distance);
    }
    return min_distance;
}

double
getMinimumDistanceBetweenSpheres(const std::vector<std::pair<Vec3d, double>>& spheres)
{
    double min_distance = std::numeric_limits<double>::infinity();
    for (auto it = spheres.begin(); it != spheres.end(); ++it)
    {
        const auto& [center, radius] = *it;
        min_distance =
          std::min(getSphereDistanceToSpheres(
                     center, radius, boost::make_iterator_range(it + 1, spheres.end())),
                   min_distance);
    }
    return min_distance;
}

double
getMaximumDistanceOfOneSphereToOthers(const std::vector<std::pair<Vec3d, double>>& spheres)
{
    double max_distance = -std::numeric_limits<double>::infinity();
    for (auto it = spheres.begin(); it != spheres.end(); ++it)
    {
        const auto& [center, radius] = *it;
        auto is_not_sphere = [&it](const std::pair<Vec3d, double>& sphere) {
            return sphere != *it;
        };
        max_distance =
          std::max(getSphereDistanceToSpheres(
                     center, radius, spheres | boost::adaptors::filtered(is_not_sphere)),
                   max_distance);
    }
    return max_distance;
}

bool
isSphereIntersectingOtherSpheres(const Vec3d& center,
                                 double radius,
                                 const std::vector<std::pair<Vec3d, double>>& other_spheres,
                                 double tolerance = 0.01)
{

    return getSphereDistanceToSpheres(center, radius, other_spheres) < tolerance;
}

// TODO: Is copying the distribution the best ?
template<class RadiusDistribution, class CoordinateDistribution>
std::enable_if_t<
  std::is_same_v<typename RadiusDistribution::result_type,
                 double> && std::is_same_v<typename CoordinateDistribution::result_type, double>,
  std::vector<std::pair<Vec3d, double>>>
getRandomSpheres(RadiusDistribution radius_distribution,
                 CoordinateDistribution coordinate_distribution,
                 std::size_t number_spheres,
                 std::size_t max_number_of_tries,
                 bool non_intersecting = false)
{
    std::default_random_engine engine(std::random_device{}());

    std::vector<std::pair<Vec3d, double>> spheres;
    for (std::size_t i = 0; i < number_spheres; ++i)
    {
        for (std::size_t j = 0; j < max_number_of_tries; ++j)
        {
            Vec3d center = Vec3d(coordinate_distribution(engine),
                                 coordinate_distribution(engine),
                                 coordinate_distribution(engine));
            double radius = radius_distribution(engine);

            if (non_intersecting && isSphereIntersectingOtherSpheres(center, radius, spheres))
            {
                continue;
            }

            spheres.push_back(std::pair<Vec3d, double>(center, radius));
            break;
        }
    }

    return spheres;
}

// TODO: Is copying the distribution the best ?
template<class RadiusDistribution, class CoordinateDistribution>
std::enable_if_t<
  std::is_same_v<typename RadiusDistribution::result_type,
                 double> && std::is_same_v<typename CoordinateDistribution::result_type, double>,
  std::vector<std::pair<Vec3d, double>>>
getRandomNonIntersectingSpheres(RadiusDistribution radius_distribution,
                                CoordinateDistribution coordinate_distribution,
                                std::size_t number_spheres,
                                std::size_t max_number_of_tries)
{
    return getRandomSpheres(
      radius_distribution, coordinate_distribution, number_spheres, max_number_of_tries, true);
}

#warning Factorize this into one templated function
void
convertToLosTopos(const std::vector<Vec3d>& in, std::vector<LosTopos::Vec3d>& out)
{
    out.resize(in.size());
    for (size_t i = 0; i < in.size(); i++)
    {
        out[i] = LosTopos::Vec3d(in[i][0], in[i][1], in[i][2]);
    }
}

void
convertToLosTopos(const std::vector<Vec3i>& in, std::vector<LosTopos::Vec3st>& out)
{
    out.resize(in.size());
    for (size_t i = 0; i < in.size(); i++)
    {
        out[i] = LosTopos::Vec3st(in[i][0], in[i][1], in[i][2]);
    }
}

void
convertToLosTopos(const std::vector<Vec2i>& in, std::vector<LosTopos::Vec2i>& out)
{
    out.resize(in.size());
    for (size_t i = 0; i < in.size(); i++)
    {
        out[i] = LosTopos::Vec2i(in[i][0], in[i][1]);
    }
}
}

VS3D*
Scenes::sceneSphere(Sim* sim,
                    std::vector<LosTopos::Vec3d>& vs,
                    std::vector<LosTopos::Vec3st>& fs,
                    std::vector<LosTopos::Vec2i>& ls,
                    std::vector<size_t>& cv,
                    std::vector<Vec3d>& cx)
{
    int N = Options::intValue("mesh-size-n");
    //    int M = Options::intValue("mesh-size-m");

    std::vector<Vec3d> v;
    std::vector<Vec3i> f;
    std::vector<Vec2i> l;

    double r = 1.0;
    createIcoSphere(Vec3d(0, 0, 0), r, N, v, f, l);

    vs.resize(v.size());
    fs.resize(f.size());
    ls.resize(l.size());
    for (size_t i = 0; i < v.size(); i++)
        vs[i] = LosTopos::Vec3d(v[i][0] * 0.4, v[i][1], v[i][2]);
    for (size_t i = 0; i < f.size(); i++)
        fs[i] = LosTopos::Vec3st(f[i][0], f[i][1], f[i][2]);
    for (size_t i = 0; i < l.size(); i++)
        ls[i] = LosTopos::Vec2i(l[i][0], l[i][1]);

    return new VS3D(vs, fs, ls, cv, cx);
}

VS3D*
Scenes::sceneTet(Sim* sim,
                 std::vector<LosTopos::Vec3d>& vs,
                 std::vector<LosTopos::Vec3st>& fs,
                 std::vector<LosTopos::Vec2i>& ls,
                 std::vector<size_t>& cv,
                 std::vector<Vec3d>& cx)
{
    int N = Options::intValue("mesh-size-n");

    std::vector<Vec3d> v;
    std::vector<Vec3i> f;
    std::vector<Vec2i> l;

    double r = 1.0;
    double len = sqrt(3.0);
    v.push_back(Vec3d(-1, -1, -1) * r / len);
    v.push_back(Vec3d(-1, 1, 1) * r / len);
    v.push_back(Vec3d(1, -1, 1) * r / len);
    v.push_back(Vec3d(1, 1, -1) * r / len);

    f.push_back(Vec3i(0, 2, 1));
    l.push_back(Vec2i(1, 0));
    f.push_back(Vec3i(0, 3, 2));
    l.push_back(Vec2i(1, 0));
    f.push_back(Vec3i(0, 1, 3));
    l.push_back(Vec2i(1, 0));
    f.push_back(Vec3i(1, 2, 3));
    l.push_back(Vec2i(1, 0));

    for (int i = 0; i < N; i++)
        subdivide(Vec3d(0, 0, 0), 0, v, f, l);

    vs.resize(v.size());
    fs.resize(f.size());
    ls.resize(l.size());
    for (size_t i = 0; i < v.size(); i++)
        vs[i] = LosTopos::Vec3d(v[i][0], v[i][1], v[i][2]);
    for (size_t i = 0; i < f.size(); i++)
        fs[i] = LosTopos::Vec3st(f[i][0], f[i][1], f[i][2]);
    for (size_t i = 0; i < l.size(); i++)
        ls[i] = LosTopos::Vec2i(l[i][0], l[i][1]);

    return new VS3D(vs, fs, ls, cv, cx);
}

VS3D*
Scenes::sceneCube(Sim* sim,
                  std::vector<LosTopos::Vec3d>& vs,
                  std::vector<LosTopos::Vec3st>& fs,
                  std::vector<LosTopos::Vec2i>& ls,
                  std::vector<size_t>& cv,
                  std::vector<Vec3d>& cx)
{
    int N = Options::intValue("mesh-size-n");

    std::vector<Vec3d> v;
    std::vector<Vec3i> f;
    std::vector<Vec2i> l;

    double r = 1.0;
    double len = sqrt(3.0);
    v.push_back(Vec3d(-1, -1, -1) * r / len);
    v.push_back(Vec3d(1, -1, -1) * r / len);
    v.push_back(Vec3d(1, 1, -1) * r / len);
    v.push_back(Vec3d(-1, 1, -1) * r / len);
    v.push_back(Vec3d(-1, -1, 1) * r / len);
    v.push_back(Vec3d(1, -1, 1) * r / len);
    v.push_back(Vec3d(1, 1, 1) * r / len);
    v.push_back(Vec3d(-1, 1, 1) * r / len);

    f.push_back(Vec3i(0, 3, 2));
    l.push_back(Vec2i(1, 0));
    f.push_back(Vec3i(2, 1, 0));
    l.push_back(Vec2i(1, 0));
    f.push_back(Vec3i(0, 1, 5));
    l.push_back(Vec2i(1, 0));
    f.push_back(Vec3i(5, 4, 0));
    l.push_back(Vec2i(1, 0));
    f.push_back(Vec3i(1, 2, 6));
    l.push_back(Vec2i(1, 0));
    f.push_back(Vec3i(6, 5, 1));
    l.push_back(Vec2i(1, 0));
    f.push_back(Vec3i(2, 3, 7));
    l.push_back(Vec2i(1, 0));
    f.push_back(Vec3i(7, 6, 2));
    l.push_back(Vec2i(1, 0));
    f.push_back(Vec3i(3, 0, 4));
    l.push_back(Vec2i(1, 0));
    f.push_back(Vec3i(4, 7, 3));
    l.push_back(Vec2i(1, 0));
    f.push_back(Vec3i(4, 5, 6));
    l.push_back(Vec2i(1, 0));
    f.push_back(Vec3i(6, 7, 4));
    l.push_back(Vec2i(1, 0));

    for (int i = 0; i < N; i++)
        subdivide(Vec3d(0, 0, 0), 0, v, f, l);

    vs.resize(v.size());
    fs.resize(f.size());
    ls.resize(l.size());
    for (size_t i = 0; i < v.size(); i++)
        vs[i] = LosTopos::Vec3d(v[i][0], v[i][1], v[i][2]);
    for (size_t i = 0; i < f.size(); i++)
        fs[i] = LosTopos::Vec3st(f[i][0], f[i][1], f[i][2]);
    for (size_t i = 0; i < l.size(); i++)
        ls[i] = LosTopos::Vec2i(l[i][0], l[i][1]);

    return new VS3D(vs, fs, ls, cv, cx);
}

VS3D*
Scenes::sceneSheet(Sim* sim,
                   std::vector<LosTopos::Vec3d>& vs,
                   std::vector<LosTopos::Vec3st>& fs,
                   std::vector<LosTopos::Vec2i>& ls,
                   std::vector<size_t>& cv,
                   std::vector<Vec3d>& cx)
{
    int N = Options::intValue("mesh-size-n");

    std::vector<Vec3d> v;
    std::vector<Vec3i> f;
    std::vector<Vec2i> l;

    // square sheet
    //    v.push_back(Vec3d(-1,  0, -1));
    //    v.push_back(Vec3d( 1,  0, -1));
    //    v.push_back(Vec3d(-1,  0,  1));
    //    v.push_back(Vec3d( 1,  0,  1));
    //
    //    f.push_back(Vec3i(0, 1, 3));    l.push_back(Vec2i(0, 1));
    //    f.push_back(Vec3i(3, 2, 0));    l.push_back(Vec2i(0, 1));

    // hexagonal sheet
    double h = std::sqrt(3.0) / 2;
    v.push_back(Vec3d(0, 0.2, 0));
    v.push_back(Vec3d(-1, 0.2, 0));
    v.push_back(Vec3d(-0.5, 0.2, -h));
    v.push_back(Vec3d(0.5, 0.2, -h));
    v.push_back(Vec3d(1, 0.2, 0));
    v.push_back(Vec3d(0.5, 0.2, h));
    v.push_back(Vec3d(-0.5, 0.2, h));

    f.push_back(Vec3i(0, 1, 2));
    l.push_back(Vec2i(1, 0));
    f.push_back(Vec3i(0, 2, 3));
    l.push_back(Vec2i(1, 0));
    f.push_back(Vec3i(0, 3, 4));
    l.push_back(Vec2i(1, 0));
    f.push_back(Vec3i(0, 4, 5));
    l.push_back(Vec2i(1, 0));
    f.push_back(Vec3i(0, 5, 6));
    l.push_back(Vec2i(1, 0));
    f.push_back(Vec3i(0, 6, 1));
    l.push_back(Vec2i(1, 0));

    for (int i = 0; i < N; i++)
        subdivide(Vec3d(0, 0, 0), 0, v, f, l);

    int buf[] = { 1, 26, 13, 29, 2, 30, 15, 32, 3, 33, 16, 35,
                  4, 36, 17, 38, 5, 39, 18, 42, 6, 41, 14, 27 };
    int nbuf = sizeof(buf) / sizeof(int);

    std::vector<Vec3d> vel;
    std::vector<unsigned char> fixed;
    for (int i = 0; i < nbuf; ++i)
    {
        fixed.push_back(buf[i] == 18);
        vel.push_back(Vec3d::Zero());
        cv.push_back(buf[i]);
        v[buf[i]] += Vec3d(0, -0.6 * v[buf[i]](2) * v[buf[i]](2), 0);
        cx.push_back(v[buf[i]]);
    }

    double displacement = 0.1;
    for (size_t i = 0; i < v.size(); i++)
        if ((v[i] - Vec3d(0, 0, 0)).norm() < 1e-6)
            v[i].y() += displacement;

    vs.resize(v.size());
    fs.resize(f.size());
    ls.resize(l.size());
    for (size_t i = 0; i < v.size(); i++)
        vs[i] = LosTopos::Vec3d(v[i][0], v[i][1], v[i][2]);
    for (size_t i = 0; i < f.size(); i++)
        fs[i] = LosTopos::Vec3st(f[i][0], f[i][1], f[i][2]);
    for (size_t i = 0; i < l.size(); i++)
        ls[i] = LosTopos::Vec2i(l[i][0], l[i][1]);

    return new VS3D(vs, fs, ls, cv, cx, vel, fixed);
}

VS3D*
Scenes::sceneBarrel(Sim* sim,
                    std::vector<LosTopos::Vec3d>& vs,
                    std::vector<LosTopos::Vec3st>& fs,
                    std::vector<LosTopos::Vec2i>& ls,
                    std::vector<size_t>& cv,
                    std::vector<Vec3d>& cx)
{
    int N = Options::intValue("mesh-size-n");
    int M = Options::intValue("mesh-size-m");

    std::vector<Vec3d> v;
    std::vector<Vec3i> f;
    std::vector<Vec2i> l;

    double r = 1;
    double L = 1;
    for (int i = 0; i < M; i++)
        for (int j = 0; j < N; j++)
            v.push_back(Vec3d(r * cos(j * 2 * M_PI / N),
                              ((double)i / (M - 1) - 0.5) * L,
                              -r * sin(j * 2 * M_PI / N)));

    for (int i = 0; i < M - 1; i++)
        for (int j = 0; j < N; j++)
        {
            f.push_back(Vec3i(i * N + j, i * N + (j + 1) % N, (i + 1) * N + (j + 1) % N));
            l.push_back(Vec2i(1, 0));
            f.push_back(Vec3i((i + 1) * N + (j + 1) % N, (i + 1) * N + j, i * N + j));
            l.push_back(Vec2i(1, 0));
        }

    vs.resize(v.size());
    fs.resize(f.size());
    ls.resize(l.size());
    for (size_t i = 0; i < v.size(); i++)
        vs[i] = LosTopos::Vec3d(v[i][0], v[i][1], v[i][2]);
    for (size_t i = 0; i < f.size(); i++)
        fs[i] = LosTopos::Vec3st(f[i][0], f[i][1], f[i][2]);
    for (size_t i = 0; i < l.size(); i++)
        ls[i] = LosTopos::Vec2i(l[i][0], l[i][1]);

    return new VS3D(vs, fs, ls, cv, cx);
}

VS3D*
Scenes::sceneDoubleBubble(Sim* sim,
                          std::vector<LosTopos::Vec3d>& vs,
                          std::vector<LosTopos::Vec3st>& fs,
                          std::vector<LosTopos::Vec2i>& ls,
                          std::vector<size_t>& cv,
                          std::vector<Vec3d>& cx)
{
    int N = Options::intValue("mesh-size-n");
    int M = Options::intValue("mesh-size-m");

    double r = 1.0; // radius of the two bubbles
    int noverlap = N / 8;
    double d = 2 * r * cos(M_PI * noverlap / N);

    vs.push_back(LosTopos::Vec3d(-d / 2 - r, 0, 0)); // left cap
    for (int i = 1; i < N + 1 - noverlap; i++)
        for (int j = 0; j < M; j++)
            vs.push_back(LosTopos::Vec3d(-d / 2 - r * cos(M_PI * i / N),
                                         r * sin(M_PI * i / N) * cos(2 * M_PI * j / M),
                                         r * sin(M_PI * i / N) * sin(2 * M_PI * j / M)));
    for (int i = 1; i < N - noverlap; i++)
        for (int j = 0; j < M; j++)
            vs.push_back(LosTopos::Vec3d(d / 2 + r * cos(M_PI * i / N),
                                         r * sin(M_PI * i / N) * cos(2 * M_PI * j / M),
                                         r * sin(M_PI * i / N) * sin(2 * M_PI * j / M)));
    size_t rid = vs.size();
    vs.push_back(LosTopos::Vec3d(d / 2 + r, 0, 0)); // right cap
    size_t cid = vs.size();
    vs.push_back(LosTopos::Vec3d(
      0, 0, 0)); // center vertex for triangulating the planar wall separating the two bubbles

    for (int j = 0; j < M; j++)
        fs.push_back(LosTopos::Vec3st(0, 1 + j, 1 + (j + 1) % M)),
          ls.push_back(LosTopos::Vec2i(0, 1));
    for (int i = 1; i < N - noverlap; i++)
        for (int j = 0; j < M; j++)
            fs.push_back(
              LosTopos::Vec3st(1 + (i - 1) * M + j, 1 + i * M + j, 1 + i * M + (j + 1) % M)),
              ls.push_back(LosTopos::Vec2i(0, 1)),
              fs.push_back(LosTopos::Vec3st(
                1 + i * M + (j + 1) % M, 1 + (i - 1) * M + (j + 1) % M, 1 + (i - 1) * M + j)),
              ls.push_back(LosTopos::Vec2i(0, 1));
    for (int j = 0; j < M; j++)
        fs.push_back(
          LosTopos::Vec3st(rid, 1 + (N - noverlap) * M + j, 1 + (N - noverlap) * M + (j + 1) % M)),
          ls.push_back(LosTopos::Vec2i(2, 0));
    for (int i = 1; i < N - noverlap - 1; i++)
        for (int j = 0; j < M; j++)
            fs.push_back(LosTopos::Vec3st(1 + (N - noverlap) * M + (i - 1) * M + j,
                                          1 + (N - noverlap) * M + i * M + j,
                                          1 + (N - noverlap) * M + i * M + (j + 1) % M)),
              ls.push_back(LosTopos::Vec2i(2, 0)),
              fs.push_back(LosTopos::Vec3st(1 + (N - noverlap) * M + i * M + (j + 1) % M,
                                            1 + (N - noverlap) * M + (i - 1) * M + (j + 1) % M,
                                            1 + (N - noverlap) * M + (i - 1) * M + j)),
              ls.push_back(LosTopos::Vec2i(2, 0));
    for (int j = 0; j < M; j++)
        fs.push_back(LosTopos::Vec3st(1 + (N - noverlap - 1) * 2 * M + j,
                                      1 + (N - noverlap - 1) * M + j,
                                      1 + (N - noverlap - 1) * M + (j + 1) % M)),
          ls.push_back(LosTopos::Vec2i(2, 0)),
          fs.push_back(LosTopos::Vec3st(1 + (N - noverlap - 1) * M + (j + 1) % M,
                                        1 + (N - noverlap - 1) * 2 * M + (j + 1) % M,
                                        1 + (N - noverlap - 1) * 2 * M + j)),
          ls.push_back(LosTopos::Vec2i(2, 0));
    for (int j = 0; j < M; j++)
        fs.push_back(LosTopos::Vec3st(
          cid, 1 + (N - noverlap - 1) * M + (j + 1) % M, 1 + (N - noverlap - 1) * M + j)),
          ls.push_back(LosTopos::Vec2i(2, 1));

    return new VS3D(vs, fs, ls, cv, cx);
}

VS3D*
Scenes::sceneTwoBubbles(Sim* sim,
                        std::vector<LosTopos::Vec3d>& vs,
                        std::vector<LosTopos::Vec3st>& fs,
                        std::vector<LosTopos::Vec2i>& ls,
                        std::vector<size_t>& cv,
                        std::vector<Vec3d>& cx)
{
    int N = Options::intValue("mesh-size-n");

    std::vector<Vec3d> v_all;
    std::vector<Vec3i> f_all;
    std::vector<Vec2i> l_all;

    std::vector<Vec3d> v;
    std::vector<Vec3i> f;
    std::vector<Vec2i> l;

    double r = 1.0;

    v.clear();
    f.clear();
    l.clear();
    createIcoSphere(Vec3d(-1.02, 0, 0), r, N, v, f, l, Vec2i(1, 0));
    f_all.reserve(f_all.size() + f.size());
    for (size_t i = 0; i < f.size(); i++)
        f_all.push_back(
          Vec3i(f[i][0] + v_all.size(), f[i][1] + v_all.size(), f[i][2] + v_all.size()));
    v_all.insert(v_all.end(), v.begin(), v.end());
    l_all.insert(l_all.end(), l.begin(), l.end());

    v.clear();
    f.clear();
    l.clear();
    createIcoSphere(Vec3d(1.02, 0, 0), r, N, v, f, l, Vec2i(2, 0));
    f_all.reserve(f_all.size() + f.size());
    for (size_t i = 0; i < f.size(); i++)
        f_all.push_back(
          Vec3i(f[i][0] + v_all.size(), f[i][1] + v_all.size(), f[i][2] + v_all.size()));
    v_all.insert(v_all.end(), v.begin(), v.end());
    l_all.insert(l_all.end(), l.begin(), l.end());

    vs.resize(v_all.size());
    fs.resize(f_all.size());
    ls.resize(l_all.size());
    for (size_t i = 0; i < v_all.size(); i++)
        vs[i] = LosTopos::Vec3d(v_all[i][0] * 0.7, v_all[i][1], v_all[i][2]);
    for (size_t i = 0; i < f_all.size(); i++)
        fs[i] = LosTopos::Vec3st(f_all[i][0], f_all[i][1], f_all[i][2]);
    for (size_t i = 0; i < l_all.size(); i++)
        ls[i] = LosTopos::Vec2i(l_all[i][0], l_all[i][1]);

    return new VS3D(vs, fs, ls, cv, cx);
}

VS3D*
Scenes::sceneTripleJunction(Sim* sim,
                            std::vector<LosTopos::Vec3d>& vs,
                            std::vector<LosTopos::Vec3st>& fs,
                            std::vector<LosTopos::Vec2i>& ls,
                            std::vector<size_t>& cv,
                            std::vector<Vec3d>& cx)
{
    int N = Options::intValue("mesh-size-n");

    std::vector<Vec3d> v;
    std::vector<Vec3i> f;
    std::vector<Vec2i> l;

    double len = 1.0;
    double r = 1.0;
    v.push_back(Vec3d(0, 0, 0));
    v.push_back(Vec3d(0, len, 0));
    v.push_back(Vec3d(r * 2, 0, 0));
    v.push_back(Vec3d(r * 2, len, 0));
    v.push_back(Vec3d(-r * 0.5, 0, -r * 0.866025404));
    v.push_back(Vec3d(-r * 0.5, len, -r * 0.866025404));
    v.push_back(Vec3d(-r * 0.5, 0, r * 0.866025404));
    v.push_back(Vec3d(-r * 0.5, len, r * 0.866025404));

    f.push_back(Vec3i(0, 1, 2));
    l.push_back(Vec2i(0, 1));
    f.push_back(Vec3i(3, 2, 1));
    l.push_back(Vec2i(0, 1));
    f.push_back(Vec3i(0, 1, 4));
    l.push_back(Vec2i(1, 2));
    f.push_back(Vec3i(5, 4, 1));
    l.push_back(Vec2i(1, 2));
    f.push_back(Vec3i(0, 1, 6));
    l.push_back(Vec2i(2, 0));
    f.push_back(Vec3i(7, 6, 1));
    l.push_back(Vec2i(2, 0));

    for (int i = 0; i < N; i++)
        subdivide(Vec3d(0, 0, 0), 0, v, f, l);

    vs.resize(v.size());
    fs.resize(f.size());
    ls.resize(l.size());
    for (size_t i = 0; i < v.size(); i++)
        vs[i] = LosTopos::Vec3d(v[i][0], v[i][1], v[i][2]);
    for (size_t i = 0; i < f.size(); i++)
        fs[i] = LosTopos::Vec3st(f[i][0], f[i][1], f[i][2]);
    for (size_t i = 0; i < l.size(); i++)
        ls[i] = LosTopos::Vec2i(l[i][0], l[i][1]);

    return new VS3D(vs, fs, ls, cv, cx);
}

VS3D*
Scenes::sceneFoamInit(Sim* sim,
                      std::vector<LosTopos::Vec3d>& vs,
                      std::vector<LosTopos::Vec3st>& fs,
                      std::vector<LosTopos::Vec2i>& ls,
                      std::vector<size_t>& cv,
                      std::vector<Vec3d>& cx)
{
    int M = Options::intValue("mesh-size-m"); // number of bubbles
    std::vector<std::pair<Vec3d, double>> spheres =
      getRandomSpheres(std::uniform_real_distribution(0.2, 0.5),
                       std::uniform_real_distribution(-1.0, 1.0),
                       M,
                       1000u);

    std::cout << M << " spheres requested; " << spheres.size() << " actually spheres generated."
              << std::endl;

    int N = Options::intValue("mesh-size-n"); // subdiv levels

    std::vector<Vec3d> v_all;
    std::vector<Vec3i> f_all;
    std::vector<Vec2i> l_all;

    std::vector<Vec3d> v;
    std::vector<Vec3i> f;
    std::vector<Vec2i> l;
    std::size_t sphere_interior_region = 1;
    for (const auto& [center, radius] : spheres)
    {
        v.clear();
        f.clear();
        l.clear();
        createIcoSphere(center, radius, N, v, f, l, Vec2i(sphere_interior_region, 0));
        addMesh(v, f, l, v_all, f_all, l_all);

        ++sphere_interior_region;
    }

    convertToLosTopos(v_all, vs);
    convertToLosTopos(f_all, fs);
    convertToLosTopos(l_all, ls);

    return new VS3D(vs, fs, ls, cv, cx);
}

VS3D*
Scenes::sceneFoam(Sim* sim,
                  std::vector<LosTopos::Vec3d>& vs,
                  std::vector<LosTopos::Vec3st>& fs,
                  std::vector<LosTopos::Vec2i>& ls,
                  std::vector<size_t>& cv,
                  std::vector<Vec3d>& cx)
{
    sim->m_load_directory = Options::strValue("load-dir");
    assert(sim->m_load_directory != "");

    VS3D* tmp_vs = new VS3D(std::vector<LosTopos::Vec3d>(),
                            std::vector<LosTopos::Vec3st>(),
                            std::vector<LosTopos::Vec2i>());
    MeshIO::load(*tmp_vs, sim->m_load_directory + "/mesh001000.rec");

    vs = tmp_vs->m_st->pm_positions;
    fs = tmp_vs->m_st->m_mesh.m_tris;
    ls = tmp_vs->m_st->m_mesh.m_triangle_labels;

    //    for (size_t i = 0; i < vs.size(); i++)
    //        vs[i] *= 0.5;

    return new VS3D(vs, fs, ls, cv, cx);
}

VS3D*
Scenes::sceneQuadJunction(Sim* sim,
                          std::vector<LosTopos::Vec3d>& vs,
                          std::vector<LosTopos::Vec3st>& fs,
                          std::vector<LosTopos::Vec2i>& ls,
                          std::vector<size_t>& cv,
                          std::vector<Vec3d>& cx)
{
    int N = Options::intValue("mesh-size-n");

    std::vector<Vec3d> v;
    std::vector<Vec3i> f;
    std::vector<Vec2i> l;

    double r = 1.0;
    v.push_back(Vec3d(0, 0, 0));

    v.push_back(Vec3d(-1, 0, -0.70710678) * r);
    v.push_back(Vec3d(1, 0, -0.70710678) * r);
    v.push_back(Vec3d(0, -1, 0.70710678) * r);
    v.push_back(Vec3d(0, 1, 0.70710678) * r);

    v.push_back(Vec3d(0, 0, -1.41421356) * r);
    v.push_back(Vec3d(-1, -1, 0) * r);
    v.push_back(Vec3d(-1, 1, 0) * r);
    v.push_back(Vec3d(1, -1, 0) * r);
    v.push_back(Vec3d(1, 1, 0) * r);
    v.push_back(Vec3d(0, 0, 1.41421356) * r);

    f.push_back(Vec3i(0, 1, 5));
    l.push_back(Vec2i(1, 0));
    f.push_back(Vec3i(5, 2, 0));
    l.push_back(Vec2i(1, 0));
    f.push_back(Vec3i(0, 1, 6));
    l.push_back(Vec2i(0, 2));
    f.push_back(Vec3i(6, 3, 0));
    l.push_back(Vec2i(0, 2));
    f.push_back(Vec3i(0, 1, 7));
    l.push_back(Vec2i(2, 1));
    f.push_back(Vec3i(7, 4, 0));
    l.push_back(Vec2i(2, 1));
    f.push_back(Vec3i(0, 2, 8));
    l.push_back(Vec2i(3, 0));
    f.push_back(Vec3i(8, 3, 0));
    l.push_back(Vec2i(3, 0));
    f.push_back(Vec3i(0, 2, 9));
    l.push_back(Vec2i(1, 3));
    f.push_back(Vec3i(9, 4, 0));
    l.push_back(Vec2i(1, 3));
    f.push_back(Vec3i(0, 3, 10));
    l.push_back(Vec2i(3, 2));
    f.push_back(Vec3i(10, 4, 0));
    l.push_back(Vec2i(3, 2));

    double disp = 0.7;
    v[0] += v[1] * disp;

    for (int i = 0; i < N; i++)
        subdivide(Vec3d(0, 0, 0), 0, v, f, l);

    vs.resize(v.size());
    fs.resize(f.size());
    ls.resize(l.size());
    for (size_t i = 0; i < v.size(); i++)
        vs[i] = LosTopos::Vec3d(v[i][0], v[i][1], v[i][2]);
    for (size_t i = 0; i < f.size(); i++)
        fs[i] = LosTopos::Vec3st(f[i][0], f[i][1], f[i][2]);
    for (size_t i = 0; i < l.size(); i++)
        ls[i] = LosTopos::Vec2i(l[i][0], l[i][1]);

    return new VS3D(vs, fs, ls, cv, cx);
}

VS3D*
Scenes::sceneConstrainedSphere(Sim* sim,
                               std::vector<LosTopos::Vec3d>& vs,
                               std::vector<LosTopos::Vec3st>& fs,
                               std::vector<LosTopos::Vec2i>& ls,
                               std::vector<size_t>& cv,
                               std::vector<Vec3d>& cx)
{
    int N = Options::intValue("mesh-size-n");
    //    int M = Options::intValue("mesh-size-m");

    std::vector<Vec3d> v;
    std::vector<Vec3i> f;
    std::vector<Vec2i> l;

    double r = 1.0;
    createIcoSphere(Vec3d(0, 0, 0), r, N, v, f, l);

    vs.resize(v.size());
    fs.resize(f.size());
    ls.resize(l.size());
    for (size_t i = 0; i < v.size(); i++)
        vs[i] = LosTopos::Vec3d(v[i][0] * 0.6, v[i][1], v[i][2]);
    for (size_t i = 0; i < f.size(); i++)
        fs[i] = LosTopos::Vec3st(f[i][0], f[i][1], f[i][2]);
    for (size_t i = 0; i < l.size(); i++)
        ls[i] = LosTopos::Vec2i(l[i][0], l[i][1]);

    cv.push_back(96);
    cv.push_back(101);
    cv.push_back(154);
    cv.push_back(155);
    cv.push_back(161);
    cv.push_back(160);
    cv.push_back(635);
    cv.push_back(641);
    cv.push_back(565);
    cv.push_back(566);
    cv.push_back(580);
    cv.push_back(581);
    for (size_t i = 0; i < cv.size(); i++)
        cx.push_back(vc(vs[cv[i]]));

    return new VS3D(vs, fs, ls, cv, cx);
}

VS3D*
Scenes::sceneBubbleWand(Sim* sim,
                        std::vector<LosTopos::Vec3d>& vs,
                        std::vector<LosTopos::Vec3st>& fs,
                        std::vector<LosTopos::Vec2i>& ls,
                        std::vector<size_t>& cv,
                        std::vector<Vec3d>& cx)
{
    int N = Options::intValue("mesh-size-n");

    std::vector<Vec3d> v;
    std::vector<Vec3i> f;
    std::vector<Vec2i> l;

    double r = 1.0;
    v.push_back(Vec3d(0, 0, 0));
    v.push_back(Vec3d(-1, 0, 0));
    for (int i = 0; i < N; i++)
        v.push_back(Vec3d(0, r * cos(i * 2 * M_PI / N), r * sin(i * 2 * M_PI / N)));

    for (int i = 0; i < N; i++)
    {
        f.push_back(Vec3i(0, i + 2, (i + 1) % N + 2));
        l.push_back(Vec2i(1, 0));
        f.push_back(Vec3i(1, i + 2, (i + 1) % N + 2));
        l.push_back(Vec2i(0, 1));
    }

    vs.resize(v.size());
    fs.resize(f.size());
    ls.resize(l.size());
    for (size_t i = 0; i < v.size(); i++)
        vs[i] = LosTopos::Vec3d(v[i][0], v[i][1], v[i][2]);
    for (size_t i = 0; i < f.size(); i++)
        fs[i] = LosTopos::Vec3st(f[i][0], f[i][1], f[i][2]);
    for (size_t i = 0; i < l.size(); i++)
        ls[i] = LosTopos::Vec2i(l[i][0], l[i][1]);

    cv.push_back(1);
    for (int i = 0; i < N; i++)
        cv.push_back(i + 2);
    for (size_t i = 0; i < cv.size(); i++)
        cx.push_back(vc(vs[cv[i]]));

    return new VS3D(vs, fs, ls, cv, cx);
}

VS3D*
Scenes::sceneTwoRingsPinching(Sim* sim,
                              std::vector<LosTopos::Vec3d>& vs,
                              std::vector<LosTopos::Vec3st>& fs,
                              std::vector<LosTopos::Vec2i>& ls,
                              std::vector<size_t>& cv,
                              std::vector<Vec3d>& cx)
{
    int N = Options::intValue("mesh-size-n");

    std::vector<Vec3d> v;
    std::vector<Vec3i> f;
    std::vector<Vec2i> l;

    double r = 0.5;
    double d = 0.15;
    v.push_back(Vec3d(-1 - d, 0, 0));
    v.push_back(Vec3d(1 + d, 0, 0));
    for (int i = 0; i < N; i++)
        v.push_back(Vec3d(-d, 0, 0) + Vec3d(0, cos(i * 2 * M_PI / N), sin(i * 2 * M_PI / N)) * r);
    for (int i = 0; i < N; i++)
        v.push_back(Vec3d(0, 0, 0) + Vec3d(0, cos(i * 2 * M_PI / N), sin(i * 2 * M_PI / N)) * r);
    for (int i = 0; i < N; i++)
        v.push_back(Vec3d(d, 0, 0) + Vec3d(0, cos(i * 2 * M_PI / N), sin(i * 2 * M_PI / N)) * r);

    for (int i = 0; i < N; i++)
    {
        f.push_back(Vec3i(0, i + 2, (i + 1) % N + 2));
        l.push_back(Vec2i(0, 1));
        f.push_back(Vec3i(1, i + 2 + 2 * N, (i + 1) % N + 2 + 2 * N));
        l.push_back(Vec2i(1, 0));

        f.push_back(Vec3i(i + 2, (i + 1) % N + 2, (i + 1) % N + 2 + N));
        l.push_back(Vec2i(1, 0));
        f.push_back(Vec3i((i + 1) % N + 2 + N, i + 2 + N, i + 2));
        l.push_back(Vec2i(1, 0));

        f.push_back(Vec3i(i + 2 + N, (i + 1) % N + 2 + N, (i + 1) % N + 2 + 2 * N));
        l.push_back(Vec2i(1, 0));
        f.push_back(Vec3i((i + 1) % N + 2 + 2 * N, i + 2 + 2 * N, i + 2 + N));
        l.push_back(Vec2i(1, 0));
    }

    vs.resize(v.size());
    fs.resize(f.size());
    ls.resize(l.size());
    for (size_t i = 0; i < v.size(); i++)
        vs[i] = LosTopos::Vec3d(v[i][0], v[i][1], v[i][2]);
    for (size_t i = 0; i < f.size(); i++)
        fs[i] = LosTopos::Vec3st(f[i][0], f[i][1], f[i][2]);
    for (size_t i = 0; i < l.size(); i++)
        ls[i] = LosTopos::Vec2i(l[i][0], l[i][1]);

    cv.push_back(0);
    cv.push_back(1);
    for (int i = 0; i < N; i++)
        cv.push_back(i + 2);
    for (int i = 0; i < N; i++)
        cv.push_back(i + 2 + 2 * N);
    for (size_t i = 0; i < cv.size(); i++)
        cx.push_back(vc(vs[cv[i]]));

    return new VS3D(vs, fs, ls, cv, cx);
}

VS3D*
Scenes::scenePullingFoam(Sim* sim,
                         std::vector<LosTopos::Vec3d>& vs,
                         std::vector<LosTopos::Vec3st>& fs,
                         std::vector<LosTopos::Vec2i>& ls,
                         std::vector<size_t>& cv,
                         std::vector<Vec3d>& cx)
{
    sim->m_load_directory = Options::strValue("load-dir");
    assert(sim->m_load_directory != "");

    VS3D* tmp_vs = new VS3D(std::vector<LosTopos::Vec3d>(),
                            std::vector<LosTopos::Vec3st>(),
                            std::vector<LosTopos::Vec2i>());
    MeshIO::load(*tmp_vs, sim->m_load_directory + "/mesh001000.rec");

    vs = tmp_vs->m_st->pm_positions;
    fs = tmp_vs->m_st->m_mesh.m_tris;
    ls = tmp_vs->m_st->m_mesh.m_triangle_labels;

    //    for (size_t i = 0; i < vs.size(); i++)
    //        vs[i] *= 0.5;

    for (size_t i = 0; i < vs.size(); i++)
        vs[i] = LosTopos::Vec3d(vs[i][2], vs[i][1], -vs[i][0]);

    VS3D* vs3d = new VS3D(vs, fs, ls, cv, cx);

    //    LosTopos::Vec3d left_axis(-1, -1, 0);   // for foam10
    LosTopos::Vec3d left_axis(-1, 0, 0); // for foam3
    std::vector<std::pair<int, double>> left_order;
    for (size_t i = 0; i < vs.size(); i++)
        left_order.push_back(std::pair<int, double>(i, -dot(vs[i], left_axis)));
    OrderComp comp;
    std::sort(left_order.begin(), left_order.end(), comp);

    //    LosTopos::Vec3d right_axis(1, 1, 0.5);  // for foam10
    LosTopos::Vec3d right_axis(1, 0, 0); // for foam3
    std::vector<std::pair<int, double>> right_order;
    for (size_t i = 0; i < vs.size(); i++)
        right_order.push_back(std::pair<int, double>(i, -dot(vs[i], right_axis)));
    std::sort(right_order.begin(), right_order.end(), comp);

    //    for (int i = 0; i < 30; i++)
    //        cv.push_back(left_order[i].first);
    //    for (int i = 0; i < 30; i++)
    //        cv.push_back(right_order[i].first);

    size_t leftmost = left_order[0].first;
    std::set<size_t> floodfill;
    floodfill.insert(leftmost);
    std::set<size_t> floodfill_next = floodfill;
    for (int k = 0; k < 2; k++)
    {
        floodfill = floodfill_next;
        for (std::set<size_t>::iterator ii = floodfill.begin(); ii != floodfill.end(); ii++)
        {
            size_t i = *ii;
            for (size_t j = 0; j < vs3d->mesh().m_vertex_to_edge_map[i].size(); j++)
            {
                LosTopos::Vec2st e = vs3d->mesh().m_edges[vs3d->mesh().m_vertex_to_edge_map[i][j]];
                size_t vother = (e[0] == i ? e[1] : e[0]);
                floodfill_next.insert(vother);
            }
        }
    }

    std::vector<size_t> last_ring;
    for (std::set<size_t>::iterator ii = floodfill_next.begin(); ii != floodfill_next.end(); ii++)
        if (floodfill.find(*ii) == floodfill.end())
            last_ring.push_back(*ii);

    Vec3d lr_center(0, 0, 0);
    Vec3d lr_normal(0, 0, 0);
    Mat3d lr_normal_sum = Mat3d::Zero();
    double lr_radius = 0;
    for (size_t i = 0; i < last_ring.size(); i++)
        lr_center += vs3d->pos(last_ring[i]);
    lr_center /= last_ring.size();
    for (size_t i = 0; i < last_ring.size(); i++)
    {
        for (size_t j = i + 1; j < last_ring.size(); j++)
        {
            Vec3d n =
              (vs3d->pos(last_ring[i]) - lr_center).cross(vs3d->pos(last_ring[j]) - lr_center);
            lr_normal_sum += n * n.transpose();
        }
        lr_radius += (vs3d->pos(last_ring[i]) - lr_center).norm();
    }
    lr_radius /= last_ring.size();
    Eigen::SelfAdjointEigenSolver<Mat3d> eig(lr_normal_sum);
    lr_normal = eig.eigenvectors().col(2).normalized();

    for (size_t i = 0; i < last_ring.size(); i++)
    {
        Vec3d x = vs3d->pos(last_ring[i]);
        x -= (x - lr_center).dot(lr_normal) * lr_normal;
        x = (x - lr_center).normalized() * lr_radius + lr_center;
        vs3d->m_st->pm_positions[last_ring[i]] = vs3d->m_st->pm_newpositions[last_ring[i]] = vc(x);
    }

    for (size_t i = 0; i < last_ring.size(); i++)
        cv.push_back(last_ring[i]);

    size_t rightmost = right_order[0].first;
    floodfill.clear();
    floodfill.insert(rightmost);
    floodfill_next = floodfill;
    for (int k = 0; k < 2; k++)
    {
        floodfill = floodfill_next;
        for (std::set<size_t>::iterator ii = floodfill.begin(); ii != floodfill.end(); ii++)
        {
            size_t i = *ii;
            for (size_t j = 0; j < vs3d->mesh().m_vertex_to_edge_map[i].size(); j++)
            {
                LosTopos::Vec2st e = vs3d->mesh().m_edges[vs3d->mesh().m_vertex_to_edge_map[i][j]];
                size_t vother = (e[0] == i ? e[1] : e[0]);
                floodfill_next.insert(vother);
            }
        }
    }

    last_ring.clear();
    for (std::set<size_t>::iterator ii = floodfill_next.begin(); ii != floodfill_next.end(); ii++)
        if (floodfill.find(*ii) == floodfill.end())
            last_ring.push_back(*ii);

    lr_center = Vec3d(0, 0, 0);
    lr_normal = Vec3d(0, 0, 0);
    lr_normal_sum = Mat3d::Zero();
    lr_radius = 0;
    for (size_t i = 0; i < last_ring.size(); i++)
        lr_center += vs3d->pos(last_ring[i]);
    lr_center /= last_ring.size();
    for (size_t i = 0; i < last_ring.size(); i++)
    {
        for (size_t j = i + 1; j < last_ring.size(); j++)
        {
            Vec3d n =
              (vs3d->pos(last_ring[i]) - lr_center).cross(vs3d->pos(last_ring[j]) - lr_center);
            lr_normal_sum += n * n.transpose();
        }
        lr_radius += (vs3d->pos(last_ring[i]) - lr_center).norm();
    }
    lr_radius /= last_ring.size();
    eig = Eigen::SelfAdjointEigenSolver<Mat3d>(lr_normal_sum);
    lr_normal = eig.eigenvectors().col(2).normalized();

    for (size_t i = 0; i < last_ring.size(); i++)
    {
        Vec3d x = vs3d->pos(last_ring[i]);
        x -= (x - lr_center).dot(lr_normal) * lr_normal;
        x = (x - lr_center).normalized() * lr_radius + lr_center;
        vs3d->m_st->pm_positions[last_ring[i]] = vs3d->m_st->pm_newpositions[last_ring[i]] = vc(x);
    }

    for (size_t i = 0; i < last_ring.size(); i++)
        cv.push_back(last_ring[i]);

    //    // for foam10 only
    //    constrained_vertices.push_back(414);
    //    constrained_vertices.push_back(255);
    //    constrained_vertices.erase(std::remove(constrained_vertices.begin(),
    //    constrained_vertices.end(), 1331), constrained_vertices.end());
    //    constrained_vertices.erase(std::remove(constrained_vertices.begin(),
    //    constrained_vertices.end(), 168), constrained_vertices.end());
    //    constrained_vertices.erase(std::remove(constrained_vertices.begin(),
    //    constrained_vertices.end(), 764), constrained_vertices.end());
    //    constrained_vertices.erase(std::remove(constrained_vertices.begin(),
    //    constrained_vertices.end(), 1406), constrained_vertices.end());
    //    constrained_vertices.erase(std::remove(constrained_vertices.begin(),
    //    constrained_vertices.end(), 729), constrained_vertices.end());
    //    constrained_vertices.erase(std::remove(constrained_vertices.begin(),
    //    constrained_vertices.end(), 797), constrained_vertices.end());
    //    constrained_vertices.erase(std::remove(constrained_vertices.begin(),
    //    constrained_vertices.end(), 2485), constrained_vertices.end());

    for (size_t i = 0; i < cv.size(); i++)
        cx.push_back(vs3d->pos(cv[i]));

    vs3d->m_constrained_vertices = cv;
    vs3d->m_constrained_positions = cx;

    for (size_t i = 0; i < cv.size(); i++)
        vs3d->m_st->m_masses[cv[i]] =
          LosTopos::Vec3d(1, 1, 1) * std::numeric_limits<double>::infinity();

    return vs3d;
}

VS3D*
Scenes::scenePeanutBubble(Sim* sim,
                          std::vector<LosTopos::Vec3d>& vs,
                          std::vector<LosTopos::Vec3st>& fs,
                          std::vector<LosTopos::Vec2i>& ls,
                          std::vector<size_t>& cv,
                          std::vector<Vec3d>& cx)
{
    int N = Options::intValue("mesh-size-n");
    int M = Options::intValue("mesh-size-m");

    double r = 1.0; // radius of the two bubbles
    int noverlap = N / 3;
    double d = 2 * r * cos(M_PI * noverlap / N);

    vs.push_back(LosTopos::Vec3d(-d / 2 - r, 0, 0)); // left cap
    for (int i = 1; i < N + 1 - noverlap; i++)
        for (int j = 0; j < M; j++)
            vs.push_back(LosTopos::Vec3d(-d / 2 - r * cos(M_PI * i / N),
                                         r * sin(M_PI * i / N) * cos(2 * M_PI * j / M),
                                         r * sin(M_PI * i / N) * sin(2 * M_PI * j / M)));
    double right_ratio = 0.999;
    for (int i = 1; i < N - noverlap; i++)
        for (int j = 0; j < M; j++)
            vs.push_back(LosTopos::Vec3d(d / 2, 0, 0)
                         + LosTopos::Vec3d(cos(M_PI * i / N),
                                           sin(M_PI * i / N) * cos(2 * M_PI * j / M),
                                           sin(M_PI * i / N) * sin(2 * M_PI * j / M))
                             * r * (right_ratio + (1 - right_ratio) * i / (N - noverlap)));
    size_t rid = vs.size();
    vs.push_back(LosTopos::Vec3d(d / 2 + r * right_ratio, 0, 0)); // right cap
    size_t cid = vs.size();
    vs.push_back(LosTopos::Vec3d(
      0, 0, 0)); // center vertex for triangulating the planar wall separating the two bubbles

    for (int j = 0; j < M; j++)
        fs.push_back(LosTopos::Vec3st(0, 1 + j, 1 + (j + 1) % M)),
          ls.push_back(LosTopos::Vec2i(0, 1));
    for (int i = 1; i < N - noverlap; i++)
        for (int j = 0; j < M; j++)
            fs.push_back(
              LosTopos::Vec3st(1 + (i - 1) * M + j, 1 + i * M + j, 1 + i * M + (j + 1) % M)),
              ls.push_back(LosTopos::Vec2i(0, 1)),
              fs.push_back(LosTopos::Vec3st(
                1 + i * M + (j + 1) % M, 1 + (i - 1) * M + (j + 1) % M, 1 + (i - 1) * M + j)),
              ls.push_back(LosTopos::Vec2i(0, 1));
    for (int j = 0; j < M; j++)
        fs.push_back(
          LosTopos::Vec3st(rid, 1 + (N - noverlap) * M + j, 1 + (N - noverlap) * M + (j + 1) % M)),
          ls.push_back(LosTopos::Vec2i(1, 0));
    for (int i = 1; i < N - noverlap - 1; i++)
        for (int j = 0; j < M; j++)
            fs.push_back(LosTopos::Vec3st(1 + (N - noverlap) * M + (i - 1) * M + j,
                                          1 + (N - noverlap) * M + i * M + j,
                                          1 + (N - noverlap) * M + i * M + (j + 1) % M)),
              ls.push_back(LosTopos::Vec2i(1, 0)),
              fs.push_back(LosTopos::Vec3st(1 + (N - noverlap) * M + i * M + (j + 1) % M,
                                            1 + (N - noverlap) * M + (i - 1) * M + (j + 1) % M,
                                            1 + (N - noverlap) * M + (i - 1) * M + j)),
              ls.push_back(LosTopos::Vec2i(1, 0));
    for (int j = 0; j < M; j++)
        fs.push_back(LosTopos::Vec3st(1 + (N - noverlap - 1) * 2 * M + j,
                                      1 + (N - noverlap - 1) * M + j,
                                      1 + (N - noverlap - 1) * M + (j + 1) % M)),
          ls.push_back(LosTopos::Vec2i(1, 0)),
          fs.push_back(LosTopos::Vec3st(1 + (N - noverlap - 1) * M + (j + 1) % M,
                                        1 + (N - noverlap - 1) * 2 * M + (j + 1) % M,
                                        1 + (N - noverlap - 1) * 2 * M + j)),
          ls.push_back(LosTopos::Vec2i(1, 0));

    for (size_t i = 0; i < M; i++)
        cv.push_back(1 + (N - noverlap - 1) * M + i);
    for (size_t i = 0; i < cv.size(); i++)
        cx.push_back(vc(vs[cv[i]]));

    return new VS3D(vs, fs, ls, cv, cx);
}

VS3D*
Scenes::sceneStraw(Sim* sim,
                   std::vector<LosTopos::Vec3d>& vs,
                   std::vector<LosTopos::Vec3st>& fs,
                   std::vector<LosTopos::Vec2i>& ls,
                   std::vector<size_t>& cv,
                   std::vector<Vec3d>& cx)
{
    int N = Options::intValue("mesh-size-n");

    std::vector<Vec3d> v;
    std::vector<Vec3i> f;
    std::vector<Vec2i> l;

    double r = 1.0;
    v.push_back(Vec3d(0, 0, 0));
    v.push_back(Vec3d(-1, 0, 0));
    for (int i = 0; i < N; i++)
        v.push_back(Vec3d(0, r * cos(i * 2 * M_PI / N), r * sin(i * 2 * M_PI / N)));

    for (int i = 0; i < N; i++)
    {
        f.push_back(Vec3i(0, i + 2, (i + 1) % N + 2));
        l.push_back(Vec2i(1, 0));
        f.push_back(Vec3i(1, i + 2, (i + 1) % N + 2));
        l.push_back(Vec2i(0, 1));
    }

    vs.resize(v.size());
    fs.resize(f.size());
    ls.resize(l.size());
    for (size_t i = 0; i < v.size(); i++)
        vs[i] = LosTopos::Vec3d(v[i][0], v[i][1], v[i][2]);
    for (size_t i = 0; i < f.size(); i++)
        fs[i] = LosTopos::Vec3st(f[i][0], f[i][1], f[i][2]);
    for (size_t i = 0; i < l.size(); i++)
        ls[i] = LosTopos::Vec2i(l[i][0], l[i][1]);

    cv.push_back(1);
    for (int i = 0; i < N; i++)
        cv.push_back(i + 2);
    for (size_t i = 0; i < cv.size(); i++)
        cx.push_back(vc(vs[cv[i]]));

    return new VS3D(vs, fs, ls, cv, cx);
}

VS3D*
Scenes::sceneCarousel(Sim* sim,
                      std::vector<LosTopos::Vec3d>& vs,
                      std::vector<LosTopos::Vec3st>& fs,
                      std::vector<LosTopos::Vec2i>& ls,
                      std::vector<size_t>& cv,
                      std::vector<Vec3d>& cx)
{
    int N = Options::intValue("mesh-size-n");
    int M = Options::intValue("mesh-size-m");

    // double bubble
    double r = 4.0; // radius of the two bubbles
    int noverlap = N / 3;
    double d = 2 * r * cos(M_PI * noverlap / N);

    vs.push_back(LosTopos::Vec3d(-d / 2 - r, 0, 0)); // left cap
    for (int i = 1; i < N + 1 - noverlap; i++)
        for (int j = 0; j < M; j++)
            vs.push_back(LosTopos::Vec3d(-d / 2 - r * cos(M_PI * i / N),
                                         r * sin(M_PI * i / N) * cos(2 * M_PI * j / M),
                                         r * sin(M_PI * i / N) * sin(2 * M_PI * j / M)));
    for (int i = 1; i < N - noverlap; i++)
        for (int j = 0; j < M; j++)
            vs.push_back(LosTopos::Vec3d(d / 2 + r * cos(M_PI * i / N),
                                         r * sin(M_PI * i / N) * cos(2 * M_PI * j / M),
                                         r * sin(M_PI * i / N) * sin(2 * M_PI * j / M)));
    size_t rid = vs.size();
    vs.push_back(LosTopos::Vec3d(d / 2 + r, 0, 0)); // right cap
    size_t cid = vs.size();
    vs.push_back(LosTopos::Vec3d(
      0, 0, 0)); // center vertex for triangulating the planar wall separating the two bubbles

    for (int j = 0; j < M; j++)
        fs.push_back(LosTopos::Vec3st(0, 1 + j, 1 + (j + 1) % M)),
          ls.push_back(LosTopos::Vec2i(0, 1));
    for (int i = 1; i < N - noverlap; i++)
        for (int j = 0; j < M; j++)
            fs.push_back(
              LosTopos::Vec3st(1 + (i - 1) * M + j, 1 + i * M + j, 1 + i * M + (j + 1) % M)),
              ls.push_back(LosTopos::Vec2i(0, 1)),
              fs.push_back(LosTopos::Vec3st(
                1 + i * M + (j + 1) % M, 1 + (i - 1) * M + (j + 1) % M, 1 + (i - 1) * M + j)),
              ls.push_back(LosTopos::Vec2i(0, 1));
    for (int j = 0; j < M; j++)
        fs.push_back(
          LosTopos::Vec3st(rid, 1 + (N - noverlap) * M + j, 1 + (N - noverlap) * M + (j + 1) % M)),
          ls.push_back(LosTopos::Vec2i(2, 0));
    for (int i = 1; i < N - noverlap - 1; i++)
        for (int j = 0; j < M; j++)
            fs.push_back(LosTopos::Vec3st(1 + (N - noverlap) * M + (i - 1) * M + j,
                                          1 + (N - noverlap) * M + i * M + j,
                                          1 + (N - noverlap) * M + i * M + (j + 1) % M)),
              ls.push_back(LosTopos::Vec2i(2, 0)),
              fs.push_back(LosTopos::Vec3st(1 + (N - noverlap) * M + i * M + (j + 1) % M,
                                            1 + (N - noverlap) * M + (i - 1) * M + (j + 1) % M,
                                            1 + (N - noverlap) * M + (i - 1) * M + j)),
              ls.push_back(LosTopos::Vec2i(2, 0));
    for (int j = 0; j < M; j++)
        fs.push_back(LosTopos::Vec3st(1 + (N - noverlap - 1) * 2 * M + j,
                                      1 + (N - noverlap - 1) * M + j,
                                      1 + (N - noverlap - 1) * M + (j + 1) % M)),
          ls.push_back(LosTopos::Vec2i(2, 0)),
          fs.push_back(LosTopos::Vec3st(1 + (N - noverlap - 1) * M + (j + 1) % M,
                                        1 + (N - noverlap - 1) * 2 * M + (j + 1) % M,
                                        1 + (N - noverlap - 1) * 2 * M + j)),
          ls.push_back(LosTopos::Vec2i(2, 0));
    for (int j = 0; j < M; j++)
        fs.push_back(LosTopos::Vec3st(
          cid, 1 + (N - noverlap - 1) * M + (j + 1) % M, 1 + (N - noverlap - 1) * M + j)),
          ls.push_back(LosTopos::Vec2i(2, 1));

    for (size_t i = 0; i < M; i++)
        cv.push_back(1 + (N / 6) * M + i);
    for (size_t i = 0; i < M; i++)
        cv.push_back(1 + ((N - noverlap) + N / 6) * M + i);

    // straw
    int nstart_straw = vs.size();
    N /= 2;

    Mat3d R;
    R << 0, -1, 0, 1, 0, 0, 0, 0, 1;
    Vec3d t(0, -6, 0);

    std::vector<Vec3d> v;
    std::vector<Vec3i> f;
    std::vector<Vec2i> l;

    r = 1.0;
    v.push_back(Vec3d(0, 0, 0));
    v.push_back(Vec3d(-1, 0, 0));
    for (int i = 0; i < N; i++)
        v.push_back(Vec3d(0, r * cos(i * 2 * M_PI / N), r * sin(i * 2 * M_PI / N)));

    for (size_t i = 0; i < v.size(); i++)
        v[i] = R * v[i] + t;

    for (int i = 0; i < N; i++)
    {
        f.push_back(Vec3i(0, i + 2, (i + 1) % N + 2)
                    + Vec3i(nstart_straw, nstart_straw, nstart_straw));
        l.push_back(Vec2i(3, 0));
        f.push_back(Vec3i(1, i + 2, (i + 1) % N + 2)
                    + Vec3i(nstart_straw, nstart_straw, nstart_straw));
        l.push_back(Vec2i(0, 3));
    }

    for (size_t i = 0; i < v.size(); i++)
        vs.push_back(LosTopos::Vec3d(v[i][0], v[i][1], v[i][2]));
    for (size_t i = 0; i < f.size(); i++)
        fs.push_back(LosTopos::Vec3st(f[i][0], f[i][1], f[i][2]));
    for (size_t i = 0; i < l.size(); i++)
        ls.push_back(LosTopos::Vec2i(l[i][0], l[i][1]));

    cv.push_back(nstart_straw + 1);
    for (int i = 0; i < N; i++)
        cv.push_back(nstart_straw + i + 2);

    for (size_t i = 0; i < cv.size(); i++)
        cx.push_back(vc(vs[cv[i]]));

    return new VS3D(vs, fs, ls, cv, cx);
}

VS3D*
Scenes::sceneOctahedron(Sim* sim,
                        std::vector<LosTopos::Vec3d>& vs,
                        std::vector<LosTopos::Vec3st>& fs,
                        std::vector<LosTopos::Vec2i>& ls,
                        std::vector<size_t>& cv,
                        std::vector<Vec3d>& cx)
{
    int N = Options::intValue("mesh-size-n");

    std::vector<Vec3d> v;
    std::vector<Vec3i> f;
    std::vector<Vec2i> l;

    Vec3d perturbation = 0.1
                         * Vec3d(2.0 * rand() / RAND_MAX - 1.0,
                                 2.0 * rand() / RAND_MAX - 1.0,
                                 2.0 * rand() / RAND_MAX - 1.0);

    double r = 1.0;
    v.push_back(Vec3d(0, 0, 0) * r + perturbation);

    v.push_back(Vec3d(1, 0, 0) * r);
    v.push_back(Vec3d(-1, 0, 0) * r);
    v.push_back(Vec3d(0, 1, 0) * r);
    v.push_back(Vec3d(0, -1, 0) * r);
    v.push_back(Vec3d(0, 0, 1) * r);
    v.push_back(Vec3d(0, 0, -1) * r);

    double a = 0.3;
    v.push_back(Vec3d(a, a, a) * r);
    v.push_back(Vec3d(-a, -a, a) * r);
    v.push_back(Vec3d(-a, a, -a) * r);
    v.push_back(Vec3d(a, -a, -a) * r);

    double rp = r * 0.5;
    v.push_back(Vec3d(1, 1, 1) * rp);
    v.push_back(Vec3d(1, 1, -1) * rp);
    v.push_back(Vec3d(1, -1, 1) * rp);
    v.push_back(Vec3d(1, -1, -1) * rp);
    v.push_back(Vec3d(-1, 1, 1) * rp);
    v.push_back(Vec3d(-1, 1, -1) * rp);
    v.push_back(Vec3d(-1, -1, 1) * rp);
    v.push_back(Vec3d(-1, -1, -1) * rp);

    f.push_back(Vec3i(0, 7, 8));
    l.push_back(Vec2i(1, 2));
    f.push_back(Vec3i(0, 7, 9));
    l.push_back(Vec2i(3, 1));
    f.push_back(Vec3i(0, 7, 10));
    l.push_back(Vec2i(2, 3));
    f.push_back(Vec3i(0, 8, 9));
    l.push_back(Vec2i(1, 4));
    f.push_back(Vec3i(0, 8, 10));
    l.push_back(Vec2i(4, 2));
    f.push_back(Vec3i(0, 9, 10));
    l.push_back(Vec2i(3, 4));

    f.push_back(Vec3i(5, 7, 8));
    l.push_back(Vec2i(2, 1));
    f.push_back(Vec3i(3, 7, 9));
    l.push_back(Vec2i(1, 3));
    f.push_back(Vec3i(1, 7, 10));
    l.push_back(Vec2i(3, 2));
    f.push_back(Vec3i(2, 8, 9));
    l.push_back(Vec2i(4, 1));
    f.push_back(Vec3i(4, 8, 10));
    l.push_back(Vec2i(2, 4));
    f.push_back(Vec3i(6, 9, 10));
    l.push_back(Vec2i(4, 3));

    f.push_back(Vec3i(1, 3, 7));
    l.push_back(Vec2i(3, 5));
    f.push_back(Vec3i(2, 4, 8));
    l.push_back(Vec2i(4, 6));
    f.push_back(Vec3i(2, 3, 9));
    l.push_back(Vec2i(1, 7));
    f.push_back(Vec3i(1, 4, 10));
    l.push_back(Vec2i(2, 8));

    f.push_back(Vec3i(1, 5, 7));
    l.push_back(Vec2i(5, 2));
    f.push_back(Vec3i(3, 5, 7));
    l.push_back(Vec2i(1, 5));
    f.push_back(Vec3i(2, 5, 8));
    l.push_back(Vec2i(6, 1));
    f.push_back(Vec3i(4, 5, 8));
    l.push_back(Vec2i(2, 6));
    f.push_back(Vec3i(2, 6, 9));
    l.push_back(Vec2i(7, 4));
    f.push_back(Vec3i(3, 6, 9));
    l.push_back(Vec2i(3, 7));
    f.push_back(Vec3i(1, 6, 10));
    l.push_back(Vec2i(8, 3));
    f.push_back(Vec3i(4, 6, 10));
    l.push_back(Vec2i(4, 8));

    f.push_back(Vec3i(1, 3, 11));
    l.push_back(Vec2i(5, 0));
    f.push_back(Vec3i(3, 5, 11));
    l.push_back(Vec2i(5, 0));
    f.push_back(Vec3i(5, 1, 11));
    l.push_back(Vec2i(5, 0));
    f.push_back(Vec3i(1, 3, 12));
    l.push_back(Vec2i(0, 3));
    f.push_back(Vec3i(3, 6, 12));
    l.push_back(Vec2i(0, 3));
    f.push_back(Vec3i(6, 1, 12));
    l.push_back(Vec2i(0, 3));
    f.push_back(Vec3i(1, 4, 13));
    l.push_back(Vec2i(0, 2));
    f.push_back(Vec3i(4, 5, 13));
    l.push_back(Vec2i(0, 2));
    f.push_back(Vec3i(5, 1, 13));
    l.push_back(Vec2i(0, 2));
    f.push_back(Vec3i(1, 4, 14));
    l.push_back(Vec2i(8, 0));
    f.push_back(Vec3i(4, 6, 14));
    l.push_back(Vec2i(8, 0));
    f.push_back(Vec3i(6, 1, 14));
    l.push_back(Vec2i(8, 0));
    f.push_back(Vec3i(2, 3, 15));
    l.push_back(Vec2i(0, 1));
    f.push_back(Vec3i(3, 5, 15));
    l.push_back(Vec2i(0, 1));
    f.push_back(Vec3i(5, 2, 15));
    l.push_back(Vec2i(0, 1));
    f.push_back(Vec3i(2, 3, 16));
    l.push_back(Vec2i(7, 0));
    f.push_back(Vec3i(3, 6, 16));
    l.push_back(Vec2i(7, 0));
    f.push_back(Vec3i(6, 2, 16));
    l.push_back(Vec2i(7, 0));
    f.push_back(Vec3i(2, 4, 17));
    l.push_back(Vec2i(6, 0));
    f.push_back(Vec3i(4, 5, 17));
    l.push_back(Vec2i(6, 0));
    f.push_back(Vec3i(5, 2, 17));
    l.push_back(Vec2i(6, 0));
    f.push_back(Vec3i(2, 4, 18));
    l.push_back(Vec2i(0, 4));
    f.push_back(Vec3i(4, 6, 18));
    l.push_back(Vec2i(0, 4));
    f.push_back(Vec3i(6, 2, 18));
    l.push_back(Vec2i(0, 4));

    for (int i = 0; i < N; i++)
        subdivide(Vec3d(0, 0, 0), 0, v, f, l);

    vs.resize(v.size());
    fs.resize(f.size());
    ls.resize(l.size());
    for (size_t i = 0; i < v.size(); i++)
        vs[i] = LosTopos::Vec3d(v[i][0], v[i][1], v[i][2]);
    for (size_t i = 0; i < f.size(); i++)
        fs[i] = LosTopos::Vec3st(f[i][0], f[i][1], f[i][2]);
    for (size_t i = 0; i < l.size(); i++)
        ls[i] = LosTopos::Vec2i(l[i][0], l[i][1]);

    for (int i = 1; i < 7; i++)
        cv.push_back(i);
    for (int i = 11; i < 19; i++)
        cv.push_back(i);

    for (size_t i = 0; i < cv.size(); i++)
        cx.push_back(vc(vs[cv[i]]));

    return new VS3D(vs, fs, ls, cv, cx);
}

VS3D*
Scenes::sceneBubbleLattice(Sim* sim,
                           std::vector<LosTopos::Vec3d>& vs,
                           std::vector<LosTopos::Vec3st>& fs,
                           std::vector<LosTopos::Vec2i>& ls,
                           std::vector<size_t>& cv,
                           std::vector<Vec3d>& cx)
{
    std::size_t number_subdivision = Options::intValue("mesh-size-n");
    std::size_t lattice_bubble_size = Options::intValue("mesh-size-m");
    double bubble_size = 1.;
    double distance_between_bubble = 0.1;

    std::vector<Vec3d> cube_vertices;
    std::vector<Vec3i> cube_faces;
    std::vector<Vec2i> cube_labels;

    std::vector<Vec3d> vertices;
    std::vector<Vec3i> faces;
    std::vector<Vec2i> labels;

    Vec3d position;
    for (std::size_t i = 0; i < lattice_bubble_size; ++i)
    {
        for (std::size_t j = 0; j < lattice_bubble_size; ++j)
        {
            for (std::size_t k = 0; k < lattice_bubble_size; ++k)
            {
                position = (bubble_size + distance_between_bubble) * Vec3d(i, j, k);
                createCube(position,
                           bubble_size,
                           number_subdivision,
                           cube_vertices,
                           cube_faces,
                           cube_labels,
                           Vec2i(0,
                                 i + lattice_bubble_size * j
                                   + lattice_bubble_size * lattice_bubble_size * k + 1));
                addMesh(cube_vertices, cube_faces, cube_labels, vertices, faces, labels);
            }
        }
    }

    convertToLosTopos(vertices, vs);
    convertToLosTopos(faces, fs);
    convertToLosTopos(labels, ls);

    return new VS3D(vs, fs, ls, cv, cx);
}

class MergedCubicBubblesFactory
{
  public:
    MergedCubicBubblesFactory(
      const std::vector<std::vector<std::vector<bool>>>& merged_cubic_bubble_position,
      double bubble_size)
    {
        initializeRegion(merged_cubic_bubble_position);
        initializeVertices(bubble_size);
        initializeTriangles();
        cleanVertices();
    }

    void subdivide(size_t number_subdivisions)
    {
        for (size_t i : boost::irange(0lu, number_subdivisions))
        {
            ::subdivide(Vec3d(0, 0, 0), 0, m_vertices, m_triangles, m_triangles_labels);
        }
    }

    void get(std::vector<LosTopos::Vec3d>& vertices,
             std::vector<LosTopos::Vec3st>& triangles,
             std::vector<LosTopos::Vec2i>& triangles_labels) const
    {
        convertToLosTopos(m_vertices, vertices);
        convertToLosTopos(m_triangles, triangles);
        convertToLosTopos(m_triangles_labels, triangles_labels);
    }

  private:
    void initializeRegion(
      const std::vector<std::vector<std::vector<bool>>>& merged_cubic_bubble_position)
    {
        m_i_size = 0;
        m_j_size = 0;
        m_k_size = 0;

        m_i_size = merged_cubic_bubble_position.size();
        m_region.resize(m_i_size);
        int region = 1;
        for (size_t i : boost::irange(0lu, merged_cubic_bubble_position.size()))
        {
            if (m_j_size == 0)
            {
                m_j_size = merged_cubic_bubble_position[i].size();
            }
            assert(m_j_size == merged_cubic_bubble_position[i].size());

            m_region[i].resize(m_j_size);
            for (size_t j : boost::irange(0lu, merged_cubic_bubble_position[i].size()))
            {
                if (m_k_size == 0)
                {
                    m_k_size = merged_cubic_bubble_position[i].size();
                }
                assert(m_k_size == merged_cubic_bubble_position[i].size());

                m_region[i][j].resize(m_k_size, 0);
                for (size_t k : boost::irange(0lu, merged_cubic_bubble_position[i][j].size()))
                {
                    if (merged_cubic_bubble_position[i][j][k])
                    {
                        m_region[i][j][k] = region;
                        ++region;
                    }
                }
            }
        }
    }

    void initializeVertices(double bubble_size)
    {
        for (size_t i : boost::irange(0lu, m_i_size + 1))
        {
            for (size_t j : boost::irange(0lu, m_j_size + 1))
            {
                for (size_t k : boost::irange(0lu, m_k_size + 1))
                {
                    m_vertices.push_back(bubble_size * Vec3d(i, j, k));
                }
            }
        }
    }

    void initializeTriangles()
    {
        size_t vertex_0, vertex_1, vertex_2, vertex_3, vertex_4, vertex_5, vertex_6;
        int region, region_i, region_j, region_k;
        for (std::size_t i : boost::irange(0lu, m_i_size + 1))
        {
            for (std::size_t j : boost::irange(0lu, m_j_size + 1))
            {
                for (std::size_t k : boost::irange(0lu, m_k_size + 1))
                {
                    vertex_0 = getVertexIndex(i, j, k);
                    vertex_1 = getVertexIndex(i + 1, j, k);
                    vertex_2 = getVertexIndex(i, j + 1, k);
                    vertex_3 = getVertexIndex(i + 1, j + 1, k);
                    vertex_4 = getVertexIndex(i, j, k + 1);
                    vertex_5 = getVertexIndex(i + 1, j, k + 1);
                    vertex_6 = getVertexIndex(i, j + 1, k + 1);

                    region = getRegion(i, j, k);
                    region_i = getINeighboorRegion(i, j, k);
                    region_j = getJNeighboorRegion(i, j, k);
                    region_k = getKNeighboorRegion(i, j, k);

                    if (region != region_k)
                    {
                        m_triangles.emplace_back(vertex_0, vertex_1, vertex_2);
                        m_triangles_labels.emplace_back(region, region_k);
                        m_triangles.emplace_back(vertex_1, vertex_3, vertex_2);
                        m_triangles_labels.emplace_back(region, region_k);
                    }

                    if (region != region_i)
                    {
                        m_triangles.emplace_back(vertex_0, vertex_2, vertex_4);
                        m_triangles_labels.emplace_back(region, region_i);
                        m_triangles.emplace_back(vertex_2, vertex_6, vertex_4);
                        m_triangles_labels.emplace_back(region, region_i);
                    }

                    if (region != region_j)
                    {
                        m_triangles.emplace_back(vertex_0, vertex_4, vertex_5);
                        m_triangles_labels.emplace_back(region, region_j);
                        m_triangles.emplace_back(vertex_0, vertex_5, vertex_1);
                        m_triangles_labels.emplace_back(region, region_j);
                    }
                }
            }
        }
    }

    void cleanVertices()
    {
        std::vector<int> mapping(m_vertices.size(), -1);
        std::vector<Vec3d> new_vertices;
        int new_index = 0;
        for (Vec3i& triangle : m_triangles)
        {
            for (size_t i : boost::irange((Vec3i::Index)0, triangle.size()))
            {
                if (mapping[triangle[i]] == -1)
                {
                    mapping[triangle[i]] = new_index;
                    new_vertices.push_back(m_vertices[triangle[i]]);
                    ++new_index;
                }
                triangle[i] = mapping[triangle[i]];
            }
        }
        m_vertices = new_vertices;
    }

    size_t getVertexIndex(size_t i, size_t j, size_t k) const
    {
        return i * (m_j_size + 1) * (m_k_size + 1) + j * (m_k_size + 1) + k;
    }

    int getINeighboorRegion(size_t i, size_t j, size_t k) const
    {
        if (i == 0)
        {
            return 0;
        }

        return getRegion(i - 1, j, k);
    }

    int getJNeighboorRegion(size_t i, size_t j, size_t k) const
    {
        if (j == 0)
        {
            return 0;
        }

        return getRegion(i, j - 1, k);
    }

    int getKNeighboorRegion(size_t i, size_t j, size_t k) const
    {
        if (k == 0)
        {
            return 0;
        }

        return getRegion(i, j, k - 1);
    }

    int getRegion(size_t i, size_t j, size_t k) const
    {
        if (i >= m_i_size || j >= m_j_size || k >= m_k_size)
        {
            return 0;
        }

        return m_region[i][j][k];
    }

    std::vector<std::vector<std::vector<int>>> m_region;
    std::vector<Vec3d> m_vertices;
    std::vector<Vec3i> m_triangles;
    std::vector<Vec2i> m_triangles_labels;
    size_t m_i_size, m_j_size, m_k_size;
};

VS3D*
Scenes::sceneMergedBubbleLattice(Sim* sim,
                                 std::vector<LosTopos::Vec3d>& vs,
                                 std::vector<LosTopos::Vec3st>& fs,
                                 std::vector<LosTopos::Vec2i>& ls,
                                 std::vector<size_t>& cv,
                                 std::vector<Vec3d>& cx)
{
    std::size_t number_subdivisions = Options::intValue("mesh-size-n");
    std::size_t number_bubbles = Options::intValue("mesh-size-m");
    std::size_t whole_foam_size = static_cast<size_t>(std::ceil(std::cbrt(number_bubbles)));
    double bubble_size = 1.;

    std::vector<std::vector<std::vector<bool>>> bubble_positions(
      whole_foam_size,
      std::vector<std::vector<bool>>(whole_foam_size, std::vector<bool>(whole_foam_size, false)));

    size_t bubble_index = 0;
    for (size_t i : boost::irange(0lu, whole_foam_size))
    {
        for (size_t j : boost::irange(0lu, whole_foam_size))
        {
            for (size_t k : boost::irange(0lu, whole_foam_size))
            {
                if (bubble_index < number_bubbles)
                {
                    bubble_positions[i][j][k] = true;
                    ++bubble_index;
                }
            }
        }
    }

    MergedCubicBubblesFactory factory(bubble_positions, bubble_size);
    factory.subdivide(number_subdivisions);
    factory.get(vs, fs, ls);
    return new VS3D(vs, fs, ls, cv, cx);
}

VS3D*
Scenes::sceneFlyingBubbles(Sim* sim,
                           std::vector<LosTopos::Vec3d>& vs,
                           std::vector<LosTopos::Vec3st>& fs,
                           std::vector<LosTopos::Vec2i>& ls,
                           std::vector<size_t>& cv,
                           std::vector<Vec3d>& cx)
{
    std::vector<std::pair<Vec3d, double>> spheres =
      getRandomNonIntersectingSpheres(std::uniform_real_distribution(0.2, 0.5),
                                      std::uniform_real_distribution(-1.0, 1.0),
                                      Options::intValue("mesh-size-m"),
                                      1000u);

    std::vector<Vec3d> spheres_velocity;
    spheres_velocity.reserve(spheres.size());
    for (size_t velocity_index : boost::irange(0lu, spheres.size()))
    {
        spheres_velocity.push_back(Vec3d::Random());
    }

    std::vector<Vec3d> vertices;
    std::vector<Vec3i> triangles;
    std::vector<Vec2i> triangles_labels;
    std::vector<Vec3d> velocities_directions;
    std::vector<double> velocities_magnitudes;

    std::vector<Vec3d> sphere_vertices;
    std::vector<Vec3i> sphere_triangles;
    std::vector<Vec2i> sphere_triangles_labels;
    std::size_t number_subdivisions = Options::intValue("mesh-size-n");

    for (size_t sphere_index : boost::irange(0lu, spheres.size()))
    {
        auto& [center, radius] = spheres[sphere_index];

        sphere_vertices.clear();
        sphere_triangles.clear();
        sphere_triangles_labels.clear();
        createIcoSphere(center,
                        radius,
                        number_subdivisions,
                        sphere_vertices,
                        sphere_triangles,
                        sphere_triangles_labels,
                        Vec2i(sphere_index + 1, 0));

        addMesh(sphere_vertices,
                sphere_triangles,
                sphere_triangles_labels,
                vertices,
                triangles,
                triangles_labels);

        velocities_directions.insert(velocities_directions.end(),
                                     sphere_vertices.size(),
                                     spheres_velocity[sphere_index].normalized());
        velocities_magnitudes.insert(velocities_magnitudes.end(),
                                     sphere_vertices.size(),
                                     spheres_velocity[sphere_index].norm());
    }

    convertToLosTopos(vertices, vs);
    convertToLosTopos(triangles, fs);
    convertToLosTopos(triangles_labels, ls);

    return new VS3D(vs, fs, ls, velocities_directions, velocities_magnitudes, cv, cx);
}

VS3D*
Scenes::sceneBubbleLine(Sim* sim,
                        std::vector<LosTopos::Vec3d>& vs,
                        std::vector<LosTopos::Vec3st>& fs,
                        std::vector<LosTopos::Vec2i>& ls,
                        std::vector<size_t>& cv,
                        std::vector<Vec3d>& cx)
{
    size_t number_bubbles = Options::intValue("mesh-size-m");
    double bubble_size = 1.;
    double noise_coefficient = .25;

    std::vector<Vec3d> vertices;
    for (size_t k : boost::irange(0lu, number_bubbles + 1))
    {
        for (size_t j : boost::irange(0lu, 2lu))
        {
            for (size_t i : boost::irange(0lu, 2lu))
            {
                vertices.push_back(bubble_size
                                   * (Vec3d(i, j, k) + noise_coefficient * Vec3d::Random()));
            }
        }
    }

    std::vector<Vec3i> triangles;
    std::vector<Vec2i> triangles_labels;
    size_t vertex_0, vertex_1, vertex_2, vertex_3, vertex_4, vertex_5, vertex_6, vertex_7;
    for (size_t k : boost::irange(0lu, number_bubbles))
    {
        vertex_0 = k * 2 * 2;
        vertex_1 = vertex_0 + 1;
        vertex_2 = vertex_0 + 1 * 2;
        vertex_3 = vertex_0 + 1 + 1 * 2;
        vertex_4 = vertex_0 + 2 * 2;
        vertex_5 = vertex_0 + 1 + 2 * 2;
        vertex_6 = vertex_0 + 1 * 2 + 2 * 2;
        vertex_7 = vertex_0 + 1 + 1 * 2 + 2 * 2;

        triangles.emplace_back(vertex_0, vertex_1, vertex_3);
        triangles_labels.emplace_back(k + 1, k);
        triangles.emplace_back(vertex_0, vertex_3, vertex_2);
        triangles_labels.emplace_back(k + 1, k);

        triangles.emplace_back(vertex_1, vertex_3, vertex_5);
        triangles_labels.emplace_back(0, k + 1);
        triangles.emplace_back(vertex_5, vertex_3, vertex_7);
        triangles_labels.emplace_back(0, k + 1);

        triangles.emplace_back(vertex_1, vertex_5, vertex_0);
        triangles_labels.emplace_back(0, k + 1);
        triangles.emplace_back(vertex_0, vertex_5, vertex_4);
        triangles_labels.emplace_back(0, k + 1);

        triangles.emplace_back(vertex_0, vertex_4, vertex_2);
        triangles_labels.emplace_back(0, k + 1);
        triangles.emplace_back(vertex_4, vertex_6, vertex_2);
        triangles_labels.emplace_back(0, k + 1);

        triangles.emplace_back(vertex_3, vertex_2, vertex_7);
        triangles_labels.emplace_back(0, k + 1);
        triangles.emplace_back(vertex_2, vertex_6, vertex_7);
        triangles_labels.emplace_back(0, k + 1);
    }

    vertex_0 = number_bubbles * 2 * 2;
    vertex_1 = vertex_0 + 1;
    vertex_2 = vertex_0 + 1 * 2;
    vertex_3 = vertex_0 + 1 + 1 * 2;
    triangles.emplace_back(vertex_0, vertex_3, vertex_1);
    triangles_labels.emplace_back(number_bubbles, 0);
    triangles.emplace_back(vertex_0, vertex_2, vertex_3);
    triangles_labels.emplace_back(number_bubbles, 0);

    convertToLosTopos(vertices, vs);
    convertToLosTopos(triangles, fs);
    convertToLosTopos(triangles_labels, ls);

    size_t number_subdivisions = Options::intValue("mesh-size-n");
    for (std::size_t dummy : boost::irange(0lu, number_subdivisions))
    {
        subdivide(Vec3d(0, 0, 0), 0, vertices, triangles, triangles_labels);
    }

    return new VS3D(vs, fs, ls, cv, cx);
}

VS3D*
Scenes::sceneBlowingBubble(Sim* sim,
                           std::vector<LosTopos::Vec3d>& vs,
                           std::vector<LosTopos::Vec3st>& fs,
                           std::vector<LosTopos::Vec2i>& ls,
                           std::vector<size_t>& cv,
                           std::vector<Vec3d>& cx)
{
    return Scenes::sceneSphere(sim, vs, fs, ls, cv, cx);
}

static size_t
getHexagonLatticeIndex(size_t i, size_t j, size_t size_straw_lattice)
{
    return i * 3 * (size_straw_lattice + 1) + j;
}

static size_t
getHexagonLatticeSize(size_t size_straw_lattice)
{
    return (2 * size_straw_lattice + 2) * 3 * (size_straw_lattice + 1);
}

static size_t
getTranslatedHexagonVertexIndex(size_t vertex_index, size_t i, size_t j, size_t size_straw_lattice)
{
    return vertex_index + getHexagonLatticeIndex(2 * i + (j % 2), 3 * j, size_straw_lattice);
}

static size_t
getHexagonCenterIndex(size_t i, size_t j, size_t size_straw_lattice)
{
    return getTranslatedHexagonVertexIndex(
      getHexagonLatticeIndex(1, 2, size_straw_lattice), i, j, size_straw_lattice);
}

static size_t
getHexagonBackCenterIndex(size_t i, size_t j, size_t size_straw_lattice)
{
    return getHexagonLatticeSize(size_straw_lattice) + i * size_straw_lattice + j;
}

static size_t
getHexagonRegion(size_t i, size_t j, size_t size_straw_lattice)
{
    return i * size_straw_lattice + j + 1;
}

static size_t
getHexagonVertexIndex(size_t vertex_index, size_t i, size_t j, size_t size_straw_lattice)
{
    size_t left_bottom_index;
    switch (vertex_index)
    {
        case 0:
            left_bottom_index = getHexagonLatticeIndex(1, 0, size_straw_lattice);
            break;
        case 1:
            left_bottom_index = getHexagonLatticeIndex(2, 1, size_straw_lattice);
            break;
        case 2:
            left_bottom_index = getHexagonLatticeIndex(2, 3, size_straw_lattice);
            break;
        case 3:
            left_bottom_index = getHexagonLatticeIndex(1, 4, size_straw_lattice);
            break;
        case 4:
            left_bottom_index = getHexagonLatticeIndex(0, 3, size_straw_lattice);
            break;
        case 5:
            left_bottom_index = getHexagonLatticeIndex(0, 1, size_straw_lattice);
            break;
        default:
            return (size_t)-1;
    }

    return getTranslatedHexagonVertexIndex(left_bottom_index, i, j, size_straw_lattice);
}

VS3D*
Scenes::sceneStraws(Sim* sim,
                    std::vector<LosTopos::Vec3d>& vs,
                    std::vector<LosTopos::Vec3st>& fs,
                    std::vector<LosTopos::Vec2i>& ls,
                    std::vector<size_t>& cv,
                    std::vector<Vec3d>& cx)
{
    size_t number_straws = Options::intValue("mesh-size-m");
    size_t size_straw_lattice = static_cast<size_t>(std::ceil(std::sqrt(number_straws)));
    double straw_radius = 1.0;

    std::vector<Vec3d> vertices;
    std::vector<Vec3i> triangles;
    std::vector<Vec2i> triangles_labels;

    double L = std::sqrt(3) * straw_radius / 2;
    double l = straw_radius / 2;
    for (size_t i : boost::irange(0lu, 2 * size_straw_lattice + 2))
    {
        for (size_t j : boost::irange(0lu, size_straw_lattice + 1))
        {
            vertices.push_back(i * Vec3d(0, 0, L) + j * Vec3d(0, l + straw_radius, 0));
            vertices.push_back(i * Vec3d(0, 0, L) + j * Vec3d(0, l + straw_radius, 0)
                               + Vec3d(0, l, 0));
            vertices.push_back(i * Vec3d(0, 0, L) + j * Vec3d(0, l + straw_radius, 0)
                               + Vec3d(0, l + straw_radius / 2, 0));
        }
    }

    for (size_t i : boost::irange(0lu, size_straw_lattice))
    {
        for (size_t j : boost::irange(0lu, size_straw_lattice))
        {
            vertices.push_back(vertices[getHexagonCenterIndex(i, j, size_straw_lattice)]
                               + Vec3d(-straw_radius, 0, 0));
        }
    }

    std::vector<std::vector<std::pair<size_t, size_t>>> vertex_to_triangle_map(vertices.size());
    size_t bubble_index = 0;
    for (size_t i : boost::irange(0lu, size_straw_lattice))
    {
        for (size_t j : boost::irange(0lu, size_straw_lattice))
        {
            if (bubble_index >= number_straws)
            {
                continue;
            }

            for (size_t vertex_index : boost::irange(0lu, 6lu))
            {
                triangles.emplace_back(
                  getHexagonCenterIndex(i, j, size_straw_lattice),
                  getHexagonVertexIndex(vertex_index, i, j, size_straw_lattice),
                  getHexagonVertexIndex((vertex_index + 1) % 6, i, j, size_straw_lattice));
                triangles_labels.emplace_back(getHexagonRegion(i, j, size_straw_lattice), 0);
                triangles.emplace_back(
                  getHexagonBackCenterIndex(i, j, size_straw_lattice),
                  getHexagonVertexIndex(vertex_index, i, j, size_straw_lattice),
                  getHexagonVertexIndex((vertex_index + 1) % 6, i, j, size_straw_lattice));
                triangles_labels.emplace_back(0, getHexagonRegion(i, j, size_straw_lattice));
            }
            ++bubble_index;
        }
    }

    for (size_t triangle_index : boost::irange(0lu, triangles.size()))
    {
        Vec3i triangle = triangles[triangle_index];
        vertex_to_triangle_map[triangle.x()].emplace_back(triangle_index, 0);
        vertex_to_triangle_map[triangle.y()].emplace_back(triangle_index, 1);
        vertex_to_triangle_map[triangle.z()].emplace_back(triangle_index, 2);
    }

    std::vector<size_t> mapping(vertices.size(), -1);
    size_t new_index = 0;
    for (size_t vertex_index : boost::irange(0lu, vertices.size()))
    {
        if (!vertex_to_triangle_map[vertex_index].empty())
        {
            mapping[vertex_index] = new_index;
            for (auto [triangle_index, index] : vertex_to_triangle_map[vertex_index])
            {
                triangles[triangle_index][index] = new_index;
            }
            ++new_index;
        }
    }

    for (int vertex_index = vertices.size() - 1; vertex_index >= 0; --vertex_index)
    {
        if (vertex_to_triangle_map[vertex_index].empty())
        {
            vertices.erase(vertices.begin() + vertex_index);
        }
    }

    bubble_index = 0;
    for (size_t i : boost::irange(0lu, size_straw_lattice))
    {
        for (size_t j : boost::irange(0lu, size_straw_lattice))
        {
            if (bubble_index >= number_straws)
            {
                continue;
            }

            cv.push_back(mapping[getHexagonBackCenterIndex(i, j, size_straw_lattice)]);
            for (size_t vertex_index : boost::irange(0lu, 6lu))
            {
                cv.push_back(
                  mapping[getHexagonVertexIndex(vertex_index, i, j, size_straw_lattice)]);
            }
            ++bubble_index;
        }
    }
    std::sort(cv.begin(), cv.end());
    auto new_end = std::unique(cv.begin(), cv.end());
    cv.erase(new_end, cv.end());

    for (size_t i = 0; i < cv.size(); i++)
        cx.push_back(vertices[cv[i]]);

    convertToLosTopos(vertices, vs);
    convertToLosTopos(triangles, fs);
    convertToLosTopos(triangles_labels, ls);

    return new VS3D(vs, fs, ls, cv, cx);
}

void
growBubbleTowardRegionZero(VS3D& vs,
                           double total_displacement_magnitude,
                           double intermediate_displacement_magnitude)
{
    size_t number_displacement =
      std::ceil(total_displacement_magnitude / intermediate_displacement_magnitude);

    std::vector<Vec3d> vertices_normals;
    for (size_t i : boost::irange(0lu, number_displacement))
    {
        vertices_normals = vs.getVerticesNormalsTowardRegion(0);
        for (size_t vertex_index : boost::irange(0lu, vs.mesh().nv()))
        {
            vs.surfTrack()->pm_newpositions[vertex_index] =
              vs.surfTrack()->pm_positions[vertex_index]
              + intermediate_displacement_magnitude * vc(vertices_normals[vertex_index]);
        }

        double actual_displacement_magnitude;
        vs.surfTrack()->integrate(intermediate_displacement_magnitude,
                                  actual_displacement_magnitude);

        vs.improveMesh(1);
    }
};

VS3D*
Scenes::sceneNewFoam(Sim* sim,
                     std::vector<LosTopos::Vec3d>& vs,
                     std::vector<LosTopos::Vec3st>& fs,
                     std::vector<LosTopos::Vec2i>& ls,
                     std::vector<size_t>& cv,
                     std::vector<Vec3d>& cx)
{
    int M = Options::intValue("mesh-size-m"); // number of bubbles
    double bubble_radius = 1.;
    double cube_size = bubble_radius * (std::cbrt(M) / 0.56);
    std::vector<std::pair<Vec3d, double>> spheres =
      getRandomNonIntersectingSpheres(std::uniform_real_distribution(bubble_radius, bubble_radius),
                                      std::uniform_real_distribution(0., cube_size),
                                      M,
                                      1000u);

    if (spheres.size() < M)
    {
        throw std::runtime_error("Not enough spheres were generated");
    }

    int N = Options::intValue("mesh-size-n"); // subdiv levels

    std::vector<Vec3d> v_all;
    std::vector<Vec3i> f_all;
    std::vector<Vec2i> l_all;

    std::vector<Vec3d> v;
    std::vector<Vec3i> f;
    std::vector<Vec2i> l;
    std::size_t sphere_interior_region = 1;
    for (const auto& [center, radius] : spheres)
    {
        v.clear();
        f.clear();
        l.clear();
        createIcoSphere(center, radius, N, v, f, l, Vec2i(sphere_interior_region, 0));
        addMesh(v, f, l, v_all, f_all, l_all);

        ++sphere_interior_region;
    }

    convertToLosTopos(v_all, vs);
    convertToLosTopos(f_all, fs);
    convertToLosTopos(l_all, ls);

    VS3D* result = new VS3D(vs, fs, ls, cv, cx);
    growBubbleTowardRegionZero(*result,
                               0.2 * bubble_radius + getMaximumDistanceOfOneSphereToOthers(spheres),
                               0.1 * bubble_radius);
    return result;
}

VS3D*
Scenes::scene2DBubbleLattice(Sim* sim,
                             std::vector<LosTopos::Vec3d>& vs,
                             std::vector<LosTopos::Vec3st>& fs,
                             std::vector<LosTopos::Vec2i>& ls,
                             std::vector<size_t>& cv,
                             std::vector<Vec3d>& cx)
{
    std::size_t number_subdivisions = Options::intValue("mesh-size-n");
    std::size_t number_bubbles = Options::intValue("mesh-size-m");
    std::size_t whole_foam_size = static_cast<size_t>(std::ceil(std::sqrt(number_bubbles)));
    double bubble_size = 1.;

    std::vector<std::vector<std::vector<bool>>> bubble_positions(
            whole_foam_size, std::vector<std::vector<bool>>(whole_foam_size, std::vector<bool>(1, false)));
    size_t bubble_index = 0;
    for (size_t i : boost::irange(0lu, whole_foam_size))
    {
        for (size_t j : boost::irange(0lu, whole_foam_size))
        {
            if (bubble_index < number_bubbles)
            {
                bubble_positions[i][j][0] = true;
                ++bubble_index;
            }
        }
    }

    MergedCubicBubblesFactory factory(bubble_positions, bubble_size);
    factory.subdivide(number_subdivisions);
    factory.get(vs, fs, ls);
    return new VS3D(vs, fs, ls, cv, cx);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  Scene-specific time stepping
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
Scenes::stepSphere(double dt, Sim* sim, VS3D* vs)
{
}

void
Scenes::stepTet(double dt, Sim* sim, VS3D* vs)
{
}

void
Scenes::stepCube(double dt, Sim* sim, VS3D* vs)
{
}

void
Scenes::stepSheet(double dt, Sim* sim, VS3D* vs)
{
    vs->stepConstrainted(dt);
}

void
Scenes::stepBarrel(double dt, Sim* sim, VS3D* vs)
{
}

void
Scenes::stepDoubleBubble(double dt, Sim* sim, VS3D* vs)
{
}

void
Scenes::stepTwoBubbles(double dt, Sim* sim, VS3D* vs)
{
}

void
Scenes::stepTripleJunction(double dt, Sim* sim, VS3D* vs)
{
}

void
Scenes::stepFoamInit(double dt, Sim* sim, VS3D* vs)
{
    static const double alpha = 0.1;

    std::vector<Vec3d> n(vs->mesh().nt(), Vec3d(0, 0, 0));
    for (size_t i = 0; i < vs->mesh().nt(); i++)
    {
        LosTopos::Vec3st t = vs->mesh().get_triangle(i);
        Vec3d x0 = vs->pos(t[0]);
        Vec3d x1 = vs->pos(t[1]);
        Vec3d x2 = vs->pos(t[2]);
        LosTopos::Vec2i l = vs->mesh().get_triangle_label(i);
        if (l[0] == 0)
            n[i] = -(x1 - x0).cross(x2 - x0).normalized();
        else if (l[1] == 0)
            n[i] = (x1 - x0).cross(x2 - x0).normalized();
    }

    std::vector<Vec3d> displacement(vs->mesh().nv(), Vec3d(0, 0, 0));
    for (size_t i = 0; i < vs->mesh().nv(); i++)
    {
        int counter = 0;
        for (size_t j = 0; j < vs->mesh().m_vertex_to_triangle_map[i].size(); j++)
        {
            LosTopos::Vec2i l =
              vs->mesh().get_triangle_label(vs->mesh().m_vertex_to_triangle_map[i][j]);
            if (l[0] == 0 || l[1] == 0)
            {
                displacement[i] += n[vs->mesh().m_vertex_to_triangle_map[i][j]];
                counter++;
            }
        }
        if (counter == 0)
            displacement[i] = Vec3d(0, 0, 0);
        else
            displacement[i] /= counter;

        vs->surfTrack()->pm_newpositions[i] =
          vs->surfTrack()->pm_positions[i] + dt * alpha * vc(displacement[i]);
    }

    double actual_dt;
    vs->surfTrack()->integrate(dt, actual_dt);
}

void
Scenes::stepFoam(double dt, Sim* sim, VS3D* vs)
{
    // burst an outer bubble every 50 seconds
    static double s_next_burst = Options::doubleValue("foam-burst-start");
    if (sim->m_time > s_next_burst)
    {
        s_next_burst += Options::doubleValue("foam-burst-interval");
        std::cout << "Bursting a random bubble now." << std::endl;

        std::set<int> burstable_regions_set;
        for (size_t i = 0; i < vs->mesh().nt(); i++)
        {
            LosTopos::Vec2i l = vs->mesh().get_triangle_label(i);
            if (l[0] == 0 || l[1] == 0)
                burstable_regions_set.insert(l[0] == 0 ? l[1] : l[0]);
        }

        std::vector<int> burstable_regions;
        burstable_regions.assign(burstable_regions_set.begin(), burstable_regions_set.end());
        std::cout << "Eligible regions: ";
        for (size_t i = 0; i < burstable_regions.size(); i++)
            std::cout << burstable_regions[i] << " ";
        std::cout << std::endl;

        int region_to_burst = burstable_regions[rand() % burstable_regions.size()];
        std::cout << "Region to be bursted: " << region_to_burst << std::endl;

        // first, find all the triple junctions that would become manifold curves because of the
        // deletion of the bursted region. because
        //  the other two wings around those triple junctions which were originally unrelated (i.e.
        //  having different arbitrary constants) would come together and form one continuous
        //  potential field, they need to be brought to terms with each other.
        // specifically, bursting region A will remove faces with label (0, A), and for each region
        // B that A is adjacent to, faces with labels
        //  (0, B) and (A, B) will be connected into one manifold patch.
        int A = region_to_burst;
        std::set<int> Bs;
        for (size_t i = 0; i < vs->mesh().nt(); i++)
        {
            LosTopos::Vec2i l = vs->mesh().get_triangle_label(i);
            if (l[0] == region_to_burst || l[1] == region_to_burst)
                Bs.insert(l[0] == region_to_burst ? l[1] : l[0]);
        }

        for (std::set<int>::iterator Bi = Bs.begin(); Bi != Bs.end(); Bi++)
        {
            int B = *Bi;

            double mean_Gamma_0B = 0; // mean_Gamma_0B stores the mean of Gammas for region pair (0,
                                      // B) on the triple junction between (0, B) and (A, B)
            double mean_Gamma_AB = 0; // mean_Gamma_0A stores the mean of Gammas for region pair (A,
                                      // B) on the triple junction between (0, B) and (A, B)
            int mean_counter = 0;
            for (size_t i = 0; i < vs->mesh().nv(); i++)
            {
                bool incident_to_0B = false;
                bool incident_to_AB = false;
                for (size_t j = 0; j < vs->mesh().m_vertex_to_triangle_map[i].size(); j++)
                {
                    LosTopos::Vec2i l =
                      vs->mesh().get_triangle_label(vs->mesh().m_vertex_to_triangle_map[i][j]);
                    if ((l[0] == 0 && l[1] == B) || (l[1] == 0 && l[0] == B))
                        incident_to_0B = true;
                    if ((l[0] == A && l[1] == B) || (l[1] == A && l[0] == B))
                        incident_to_AB = true;
                }

                if (incident_to_0B && incident_to_AB)
                {
                    mean_Gamma_0B += (*vs->m_Gamma)[i].get(0, B);
                    mean_Gamma_AB += (*vs->m_Gamma)[i].get(A, B);
                    mean_counter++;
                }
            }

            if (mean_counter != 0)
            {
                mean_Gamma_0B /= mean_counter;
                mean_Gamma_AB /= mean_counter;
            }

            double diff = mean_Gamma_0B - mean_Gamma_AB;

            for (size_t i = 0; i < vs->mesh().nv(); i++)
            {
                bool incident_to_0B = false;
                bool incident_to_AB = false;
                for (size_t j = 0; j < vs->mesh().m_vertex_to_triangle_map[i].size(); j++)
                {
                    LosTopos::Vec2i l =
                      vs->mesh().get_triangle_label(vs->mesh().m_vertex_to_triangle_map[i][j]);
                    if ((l[0] == 0 && l[1] == B) || (l[1] == 0 && l[0] == B))
                        incident_to_0B = true;
                    if ((l[0] == A && l[1] == B) || (l[1] == A && l[0] == B))
                        incident_to_AB = true;
                }

                if (incident_to_0B && incident_to_AB)
                {
                    (*vs->m_Gamma)[i].set(
                      0, B, ((*vs->m_Gamma)[i].get(0, B) + (*vs->m_Gamma)[i].get(A, B)) / 2);
                    (*vs->m_Gamma)[i].set(A, B, 0);
                }
                else if (incident_to_0B)
                {
                    (*vs->m_Gamma)[i].set(0, B, (*vs->m_Gamma)[i].get(0, B) - diff / 2);
                }
                else if (incident_to_AB)
                {
                    (*vs->m_Gamma)[i].set(0, B, (*vs->m_Gamma)[i].get(A, B) + diff / 2);
                    (*vs->m_Gamma)[i].set(A, B, 0);
                }
            }
        }

        for (size_t i = 0; i < vs->mesh().nt(); i++)
        {
            LosTopos::Vec2i l = vs->mesh().get_triangle_label(i);
            if (l[0] == region_to_burst || l[1] == region_to_burst)
            {
                int lother = (l[0] == region_to_burst ? l[1] : l[0]);
                if (lother == 0)
                    vs->surfTrack()->remove_triangle(i);
                else
                    vs->mesh().set_triangle_label(i,
                                                  (l[0] == region_to_burst
                                                     ? LosTopos::Vec2i(0, l[1])
                                                     : LosTopos::Vec2i(l[0], 0)));
            }
        }

        vs->surfTrack()->defrag_mesh_from_scratch(vs->m_constrained_vertices);

        std::cout << "Bubble bursted." << std::endl;
    }
}

void
Scenes::stepQuadJunction(double dt, Sim* sim, VS3D* vs)
{
}

void
Scenes::stepConstrainedSphere(double dt, Sim* sim, VS3D* vs)
{
}

void
Scenes::stepBubbleWand(double dt, Sim* sim, VS3D* vs)
{
    for (size_t i = 0; i < vs->m_constrained_vertices.size(); i++)
        vs->m_constrained_positions[i] += Vec3d(-1, 0, 0) * dt;
}

void
Scenes::stepTwoRingsPinching(double dt, Sim* sim, VS3D* vs)
{
    if (true)
        for (size_t i = 0; i < vs->m_constrained_vertices.size(); i++)
            vs->m_constrained_positions[i] +=
              Vec3d((vs->m_constrained_positions[i].x() > 0 ? 1 : -1), 0, 0) * 3 * dt;
}

void
Scenes::stepPullingFoam(double dt, Sim* sim, VS3D* vs)
{
    Vec3d center0(0, 0, 0);
    int counter0 = 0;
    Vec3d center1(0, 0, 0);
    int counter1 = 0;

    for (size_t i = 0; i < vs->m_constrained_positions.size(); i++)
    {
        if (vs->m_constrained_positions[i].x() > 0)
        {
            center0 += vs->m_constrained_positions[i];
            counter0++;
        }
        else
        {
            center1 += vs->m_constrained_positions[i];
            counter1++;
        }
    }

    center0 /= counter0;
    center1 /= counter1;

    Vec3d dir = (center0 - center1).normalized();

    for (size_t i = 0; i < vs->m_constrained_vertices.size(); i++)
        vs->m_constrained_positions[i] +=
          (vs->m_constrained_positions[i].x() > 0 ? 1 : -1) * dir * 0.02 * dt;
}

void
Scenes::stepPeanutBubble(double dt, Sim* sim, VS3D* vs)
{
}

void
Scenes::stepStraw(double dt, Sim* sim, VS3D* vs)
{
    for (size_t i = 0; i < vs->m_constrained_vertices.size(); i++)
    {
        size_t cv = vs->m_constrained_vertices[i];
        (*vs->m_Gamma)[cv].set(0, 1, (*vs->m_Gamma)[cv].get(0, 1) - 10 * dt);
    }
}

void
Scenes::stepCarousel(double dt, Sim* sim, VS3D* vs)
{
    if (sim->m_time < 10)
    {
        for (size_t i = 0; i < vs->m_constrained_vertices.size(); i++)
        {
            size_t cv = vs->m_constrained_vertices[i];
            if (vs->m_constrained_positions[i][1] < -4)
                (*vs->m_Gamma)[cv].set(0, 3, (*vs->m_Gamma)[cv].get(0, 3) - 10 * dt);
        }
    }
    else
    {
        for (size_t i = 0; i < vs->m_constrained_vertices.size(); i++)
        {
            size_t cv = vs->m_constrained_vertices[i];
            if (vs->m_constrained_positions[i][1] < -4)
                vs->m_constrained_positions[i][1] -= 0.2 * dt;
        }
    }
}

void
Scenes::stepOctahedron(double dt, Sim* sim, VS3D* vs)
{
}

void
Scenes::stepBubbleLattice(double dt, Sim* sim, VS3D* vs)
{
}

void
Scenes::stepMergedBubbleLattice(double dt, Sim* sim, VS3D* vs)
{
}

void
Scenes::stepFlyingBubbles(double dt, Sim* sim, VS3D* vs)
{
}

void
Scenes::stepBubbleLine(double dt, Sim* sim, VS3D* vs)
{
}

void
Scenes::stepBlowingBubble(double dt, Sim* sim, VS3D* vs)
{

    std::vector<Vec3d> velocities = { Vec3d(0, 0, 0.05) };
    std::vector<Vec3d> positions = { Vec3d(0, 0, -1.5) };
    // double position_lattice_cell_size = 1. / 5.;
    // for (size_t i : boost::irange(0, 5))
    //{
    //    for (size_t j : boost::irange(0, 5))
    //    {
    //        positions.push_back(Vec3d(-1, -1, -1.5) + position_lattice_cell_size * Vec3d(i, j,
    //        0)); velocities.push_back(Vec3d(0, 0, 2.0));
    //    }
    //}

    vs->projectAirVelocity(velocities, positions);
}

void
Scenes::stepStraws(double dt, Sim* sim, VS3D* vs)
{
    for (size_t i = 0; i < vs->m_constrained_vertices.size(); i++)
        vs->m_constrained_positions[i] += Vec3d(-1, 0, 0) * dt;
}

void
Scenes::stepNewFoam(double dt, Sim* sim, VS3D* vs)
{
}

void
Scenes::step2DBubbleLattice(double dt, Sim* sim, VS3D* vs)
{
}

