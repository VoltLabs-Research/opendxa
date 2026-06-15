#pragma once

#include <volt/analysis/structure_analysis.h>
#include <volt/analysis/structure_analysis_context.h>
#include <volt/core/lammps_parser.h>
#include <volt/core/particle_property.h>
#include <volt/core/simulation_cell.h>
#include <volt/math/matrix3.h>

#include <memory>
#include <string>

namespace Volt{

struct FullCrystalContextParams{
    std::string referenceFile;
    int cationSpecies = 1;
    double bondCutoff = 0.0;
    std::string topologyName = "crystal";
};

struct FullCrystalContextResult{
    bool ok = false;
    std::string message;

    int atomCount = 0;
    int cationAtoms = 0;
    double grainSnapResidual = 0;

    int minNeighborCount = 0;
    int maxNeighborCount = 0;
    long atomsWithZeroNeighbors = 0;
    long edgesTotal = 0;
    long edgesWithIdealVector = 0;
    long atomsBySpeciesUnassigned = 0;

    double bulkLoopMedianResidual = -1;
    int bulkLoopsSampled = 0;

    double selectedCutoff = 0;
    double ambiguityRho = 0;

    double metricRescaleX = 1.0, metricRescaleY = 1.0, metricRescaleZ = 1.0;

    std::shared_ptr<ParticleProperty> structureTypesStorage;
};

FullCrystalContextResult buildFullCrystalContext(
    const LammpsParser::Frame& frame,
    StructureContext& context,
    StructureAnalysis& analysis,
    const FullCrystalContextParams& params
);

}