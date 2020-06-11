#include "BiotSavart.h"
#include "SimOptions.h"
#include <boost/range/algorithm/generate.hpp>
#include <boost/range/algorithm_ext/erase.hpp>
#include <boost/range/algorithm/sort.hpp>
#include <boost/range/algorithm/unique.hpp>
#include <random>
#include <boost/range/irange.hpp>
#include "fmmtl/fmmtl/util/Clock.hpp"


double
getRelativeDistance(const VecXd& a, const VecXd& b)
{
    return (a - b).cwiseAbs().maxCoeff() / (a.cwiseAbs().maxCoeff() + b.cwiseAbs().maxCoeff());
}

bool compare(const Vec3d& a, const Vec3d& b)
{
    for (size_t i : boost::irange(0lu, 3lu))
    {
        if (a[i] < b[i])
        {
            return true;
        }
        if (a[i] > b[i])
        {
            return false;
        }
    }
    return false;
}

int
main()
{
    std::random_device random_device;
    std::default_random_engine random_engine(random_device());

    std::uniform_int_distribution<size_t> number_sources_distribution(1lu, 1000lu);
    std::uniform_int_distribution<size_t> number_targets_distribution(1lu, 1000lu);
    std::uniform_real_distribution<double> coordinate_distribution(-1.,1.);
    std::uniform_real_distribution<double> charges_component_distribution(-1.,1.);

    auto generate_coordinate =
      [&coordinate_distribution, &random_engine]() {
          return coordinate_distribution(random_engine);
      };

    auto generate_position = [&generate_coordinate]() {
        return Vec3d(generate_coordinate(), generate_coordinate(), generate_coordinate());
    };

    auto generate_charge_component =
      [&charges_component_distribution, &random_engine]() {
          return charges_component_distribution(random_engine);
      };

    auto generate_charges = [&generate_charge_component]() {
        return Vec3d(
          generate_charge_component(), generate_charge_component(), generate_charge_component());
    };

    for (size_t number_sources : boost::irange(1lu, 20lu))
    {
        std::vector<Vec3d> sources(50000);
        boost::generate(sources, generate_position);
        boost::sort(sources, compare);
        boost::erase(sources, boost::unique<boost::return_found_end>(sources));

        std::vector<Vec3d> charges(sources.size());
        boost::generate(charges, generate_charges);

        std::vector<Vec3d> targets(static_cast<size_t>(std::pow(2, number_sources)));
        boost::generate(targets, generate_position);
        
        Clock computation_time;
        VecXd fast_winding_result = BiotSavart_fast_winding_number(sources, targets, charges, 0.01);
        std::cout << "ComputionTime "  << computation_time.seconds() << std::endl;
    }


    //double delta = 0.01;

    //VecXd naive_result = BiotSavart_naive(sources, targets, charges, delta);
    //VecXd fmmtl_result = BiotSavart_fmmtl(sources, targets, charges, delta);

    //Options::addDoubleOption("winding-beta", 4.);
    //Options::addIntegerOption("winding-expansion-order", 2);

    //std::cout << "Naive -- Fmmtl: " << getRelativeDistance(naive_result, fmmtl_result) << std::endl;
    //std::cout << "Fmmtl -- Fast Winding: " << getRelativeDistance(fast_winding_result, fmmtl_result)
    //          << std::endl;
    //std::cout << "Fast Winding -- Naive: " << getRelativeDistance(naive_result, fast_winding_result)
    //          << std::endl;

    return 0;
}
