#ifndef SOMA_CLUSTERING_H_
#define SOMA_CLUSTERING_H_

#include <vector>

#include "biodynamo.h"
#include "math_util.h"
#include "matrix.h"
#include "my_cell.h"
#include "soma_clustering_biology_modules.h"
#include "validation_criterion.h"

namespace bdm {

// ----------------------------------------------------------------------------
// This model examplifies the use of extracellur diffusion and shows
// how to extend the default "Cell". In step 0 one can see how an extra
// data member is added and can be accessed throughout the simulation with
// its Get and Set methods. N cells are randomly positioned in space, of which
// half are of type 1 and half of type -1.
//
// Each type secretes a different substance. Cells move towards the gradient of
// their own substance, which results in clusters being formed of cells of the
// same type.
// -----------------------------------------------------------------------------

// Define compile time parameter
template <typename Backend>
struct CompileTimeParam : public DefaultCompileTimeParam<Backend> {
  using BiologyModules = Variant<Chemotaxis, SubstanceSecretion>;
  using AtomicTypes = VariadicTypedef<MyCell>;
};

inline int Simulate(int argc, const char** argv) {
  InitializeBioDynamo(argc, argv);

  // Define initial model
  // Create an artificial bounds for the simulation space
  Param::bound_space_ = true;
  Param::min_bound_ = 0;
  Param::max_bound_ = 250;
  Param::run_mechanical_interactions_ = false;
  int num_cells = 20000;

  gTRandom.SetSeed(4357);

  // Construct num_cells/2 cells of type 1
  auto construct_0 = [](const std::array<double, 3>& position) {
    MyCell cell(position);
    cell.SetDiameter(10);
    cell.SetCellType(1);
    cell.AddBiologyModule(SubstanceSecretion());
    cell.AddBiologyModule(Chemotaxis());
    return cell;
  };
  ModelInitializer::CreateCellsRandom(Param::min_bound_, Param::max_bound_,
                                      num_cells / 2, construct_0);

  // Construct num_cells/2 cells of type -1
  auto construct_1 = [](const std::array<double, 3>& position) {
    MyCell cell(position);
    cell.SetDiameter(10);
    cell.SetCellType(-1);
    cell.AddBiologyModule(SubstanceSecretion());
    cell.AddBiologyModule(Chemotaxis());
    return cell;
  };
  ModelInitializer::CreateCellsRandom(Param::min_bound_, Param::max_bound_,
                                      num_cells / 2, construct_1);

  // Define the substances that cells may secrete
  // Order: substance_name, diffusion_coefficient, decay_constant, resolution
  ModelInitializer::DefineSubstance(kSubstance_0, "Substance_0", 0.5, 0.1, 1);
  ModelInitializer::DefineSubstance(kSubstance_1, "Substance_1", 0.5, 0.1, 1);

  // Run simulation for N timesteps
  Scheduler<> scheduler;
  scheduler.Simulate(1000);

  double spatial_range = 5;
  auto crit = GetCriterion(spatial_range, num_cells / 8);
  return !crit;
}

}  // namespace bdm

#endif  // SOMA_CLUSTERING_H_
