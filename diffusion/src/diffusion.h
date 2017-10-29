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
};

inline int Simulate(int argc, const char** argv) {
  InitializeBioDynamo(argc, argv);

  Param::backup_interval_ = 1;
  // Define initial model - in this example: two cells
  auto construct = [](const std::array<double, 3>& position) {
    Cell cell(position);
    cell.SetDiameter(30);
    cell.SetAdherence(0.4);
    cell.SetMass(1.0);
    cell.AddBiologyModule(Chemotaxis());
    // Let only one cell be responsible for the artificial substance secretion
    if (position[0] == 0 && position[1] == 0 && position[2] == 0) {
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
  ModelInitializer::CreateCells(positions, construct);

  // Define the substances that cells may secrete
  ModelInitializer::DefineSubstance(kKalium, "Kalium", 0.4, 0, 5);

  // Run simulation for N timesteps
  Param::live_visualization_ = true;
  Scheduler<> scheduler;
  scheduler.Simulate(3500);
  return 0;
}

}  // namespace bdm

#endif  // DEMO_DIFFUSION_MODULE_H_
