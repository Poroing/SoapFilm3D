#ifndef __SoapFilm3D__BiotSavart__
#define __SoapFilm3D__BiotSavart__

#include "VS3D.h"
#include "eigenheaders.h"

#include <vector>

VecXd
BiotSavart(VS3D& vs, const VecXd& dx);

VecXd
BiotSavart_naive(const std::vector<Vec3d>& sources,
                 const std::vector<Vec3d>& targets,
                 const std::vector<Vec3d>& charges,
                 double delta);

#ifndef WIN32
VecXd
BiotSavart_fmmtl(const std::vector<Vec3d>& sources,
                 const std::vector<Vec3d>& targets,
                 const std::vector<Vec3d>& charges,
                 double delta);
#endif

VecXd
BiotSavart_fast_winding_number(const std::vector<Vec3d>& sources,
                 const std::vector<Vec3d>& targets,
                 const std::vector<Vec3d>& charges,
                 double delta);

#endif //__SoapFilm3D__BiotSavart__
