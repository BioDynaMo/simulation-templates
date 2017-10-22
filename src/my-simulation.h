#ifndef MY_SIMULATION_H_
#define MY_SIMULATION_H_

#include "biodynamo.h"

namespace bdm {

// Define compile time parameter
template <typename Backend>
struct CompileTimeParam : public DefaultCompileTimeParam<Backend> {
  using BiologyModules = Variant<NullBiologyModule>;
};

inline int Simulate(int argc, const char** argv) {
  InitializeBioDynamo(argc, argv);

  // Define initial model - in this example: single cell at origin
  Cell cell({0, 0, 0});
  cell.SetDiameter(30);
  ResourceManager<>::Get()->push_back(cell);

  // Run simulation for one timestep
  Scheduler<> scheduler;
  scheduler.Simulate(1);
  return 0;
}

} // namespace bdm

#endif // MY_SIMULATION_H_