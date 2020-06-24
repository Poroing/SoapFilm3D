#include "BiotSavart.h"
#include "SimOptions.h"
#include "fmmtl/fmmtl/util/Clock.hpp"
#include <boost/range/algorithm/generate.hpp>
#include <boost/range/algorithm/sort.hpp>
#include <boost/range/algorithm/unique.hpp>
#include <boost/range/algorithm_ext/erase.hpp>
#include <boost/range/irange.hpp>
#include <random>

double
getRelativeDistance(const VecXd& a, const VecXd& b)
{
    return (a - b).cwiseAbs().maxCoeff() / (a.cwiseAbs().maxCoeff() + b.cwiseAbs().maxCoeff());
}

int
main(int argc, char** argv)
{
    std::random_device random_device;
    std::default_random_engine random_engine(random_device());

    std::uniform_real_distribution<double> coordinate_distribution(-1., 1.);
    std::uniform_real_distribution<double> charges_component_distribution(-1., 1.);

    auto generate_coordinate = [&coordinate_distribution, &random_engine]() {
        return coordinate_distribution(random_engine);
    };

    auto generate_position = [&generate_coordinate]() {
        return Vec3d(generate_coordinate(), generate_coordinate(), generate_coordinate());
    };

    auto generate_charge_component = [&charges_component_distribution, &random_engine]() {
        return charges_component_distribution(random_engine);
    };

    auto generate_charges = [&generate_charge_component]() {
        return Vec3d(
          generate_charge_component(), generate_charge_component(), generate_charge_component());
    };

    auto generate_source_position_and_charge = [&generate_position, &generate_charges]() {
        return std::make_pair(generate_position(), generate_charges());
    };

    double delta = std::atof(argv[1]);
    size_t fmmtl_expansion_order = std::atoi(argv[2]);
    double fmmtl_theta = std::atof(argv[3]);
    size_t fmmtl_minimum_cell_size = std::atoi(argv[4]);
    size_t winding_expansion_order = std::atoi(argv[5]);
    double winding_beta = std::atof(argv[6]);
    size_t number_sources = std::atoi(argv[7]);

    std::vector<std::pair<Vec3d, Vec3d>> sources_and_charges(number_sources);
    std::vector<Vec3d> targets(1000lu);

    Options::addDoubleOption("winding-beta", winding_beta);
    Options::addIntegerOption("winding-expansion-order", winding_expansion_order);
    Options::addDoubleOption("fmmtl-theta", fmmtl_theta);
    Options::addIntegerOption("fmmtl-expansion-order", fmmtl_expansion_order);
    Options::addIntegerOption("fmmtl-minimum-cell-size", fmmtl_minimum_cell_size);

    std::cout << "Fmmtl" << " " << "BarnesHut" << std::endl;
    for (size_t number_try : boost::irange(0lu, 1000lu))
    {
        boost::generate(sources_and_charges, generate_source_position_and_charge);
        erase_duplicate_sources(sources_and_charges);
        boost::generate(targets, generate_position);
        VecXd fast_winding_result =
          BiotSavart_fast_winding_number(sources_and_charges, targets, delta);
        VecXd naive_result = BiotSavart_naive(sources_and_charges, targets, delta);
        VecXd fmmtl_result = BiotSavart_fmmtl(sources_and_charges, targets, delta);
        std::cout << getRelativeDistance(naive_result, fmmtl_result) << " "
                  << getRelativeDistance(naive_result, fast_winding_result) << std::endl;
    }

    return 0;
}
