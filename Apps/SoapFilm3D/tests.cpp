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

    double delta = 0.;
    std::vector<std::pair<Vec3d, Vec3d>> sources_and_charges(
      static_cast<size_t>(std::pow(2, std::atoi(argv[1]))));
    std::vector<Vec3d> targets(1000lu);
    for (size_t number_try : boost::irange(0lu, 100lu))
    {
        boost::generate(sources_and_charges, generate_source_position_and_charge);
        erase_duplicate_sources(sources_and_charges);
        boost::generate(targets, generate_position);

        Options::addDoubleOption("winding-beta", 4.);
        Options::addIntegerOption("winding-expansion-order", 2);
        VecXd fast_winding_result =
          BiotSavart_fast_winding_number(sources_and_charges, targets, delta);
        VecXd naive_result = BiotSavart_naive(sources_and_charges, targets, delta);
        Options::addDoubleOption("fmmtl-theta", 0.5);
        Options::addIntegerOption("fmmtl-expansion-order", 20);
        Options::addIntegerOption("fmmtl-minimum-cell-size", 1);
        VecXd fmmtl_result = BiotSavart_fmmtl(sources_and_charges, targets, delta);
        std::cout << getRelativeDistance(naive_result, fmmtl_result) << " "
                  << getRelativeDistance(naive_result, fast_winding_result) << std::endl;
    }

    // double delta = 0.01;

    // std::cout << "Naive -- Fmmtl: " << getRelativeDistance(naive_result, fmmtl_result) <<
    // std::endl; std::cout << "Fmmtl -- Fast Winding: " << getRelativeDistance(fast_winding_result,
    // fmmtl_result)
    //          << std::endl;
    // std::cout << "Fast Winding -- Naive: " << getRelativeDistance(naive_result,
    // fast_winding_result)
    //          << std::endl;

    return 0;
}
