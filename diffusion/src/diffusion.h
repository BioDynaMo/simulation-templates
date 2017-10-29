#ifndef DEMO_DIFFUSION_MODULE_H_
#define DEMO_DIFFUSION_MODULE_H_

#include <vector>

#include "biodynamo.h"
#include "diffusion_biology_modules.h"

namespace bdm {

// -----------------------------------------------------------------------------
// This model creates 8 cells at each corner of a cube. A substance is
// artificially added in the middle of this cube. The cells are modeled to
// displace according to the extracellular gradient; in this case to the middle.
// -----------------------------------------------------------------------------

// Define compile time parameter
template <typename Backend>
struct CompileTimeParam : public DefaultCompileTimeParam<Backend> {
  using BiologyModules = Variant<Chemotaxis, KaliumSecretion>;
  // use default Backend and AtomicTypes
};

inline int Simulate(int argc, const char** argv) {
  // Initialize BioDynaMo
  InitializeBioDynamo(argc, argv);

  auto construct = [](const std::array<double, 3>& position) {
    Cell cell(position);
    cell.SetDiameter(30);
    cell.SetMass(1.0);
    cell.AddBiologyModule(Chemotaxis());
    std::array<double, 3> secretion_position = {{50, 50, 50}};
    if (position == secretion_position) {
      cell.AddBiologyModule(KaliumSecretion());
    }
    return cell;
  };
  std::vector<std::array<double, 3>> positions;
  positions.push_back({0, 0, 0});
  positions.push_back({100, 0, 0});
  positions.push_back({0, 100, 0});
  positions.push_back({0, 0, 100});
  positions.push_back({0, 100, 100});
  positions.push_back({100, 0, 100});
  positions.push_back({100, 100, 0});
  positions.push_back({100, 100, 100});
  // the cell responsible for secretion
  positions.push_back({50, 50, 50});
  ModelInitializer::CreateCells(positions, construct);

  // Define the substances that cells may secrete
  ModelInitializer::DefineSubstance(kKalium, "Kalium", 0.4, 0, 8);

  // Run simulation for N timesteps
  Scheduler<> scheduler;
  scheduler.Simulate(350);
  return 0;
}

}  // namespace bdm

#endif  // DEMO_DIFFUSION_MODULE_H_
