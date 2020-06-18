#ifndef __SoapFilm3D__BiotSavart__
#define __SoapFilm3D__BiotSavart__

#include "VS3D.h"
#include "eigenheaders.h"

#include <vector>

VecXd
BiotSavart(VS3D& vs, const VecXd& dx);

VecXd
BiotSavart_naive(const std::vector<std::pair<Vec3d, Vec3d>>& sources_and_charges,
                 const std::vector<Vec3d>& targets,
                 double delta);

#ifndef WIN32
VecXd
BiotSavart_fmmtl(const std::vector<std::pair<Vec3d, Vec3d>>& sources_and_charges,
                 const std::vector<Vec3d>& targets,
                 double delta);
#endif

VecXd
BiotSavart_fast_winding_number(const std::vector<std::pair<Vec3d, Vec3d>>& sources_and_charges,
                 const std::vector<Vec3d>& targets,
                 double delta);

void erase_duplicate_sources(std::vector<std::pair<Vec3d, Vec3d>>& sources_and_charges);

#endif //__SoapFilm3D__BiotSavart__
