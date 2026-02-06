#pragma once

#include <volt/core/volt.h>
#include <volt/core/frame_adapter.h>
#include <volt/core/analysis_result.h>
#include <volt/structures/crystal_structure_types.h>
#include <volt/structures/dislocation_network.h>
#include <volt/geometry/delaunay_tessellation.h>
#include <volt/analysis/elastic_mapping.h>
#include <volt/analysis/burgers_loop_builder.h>
#include <volt/geometry/interface_mesh.h>
#include <volt/math/lin_alg.h>
#include <volt/utilities/json_exporter.h>
#include <format> 
#include <vector>
#include <memory>
#include <string>

namespace Volt {

class DislocationAnalysis {
public:
    DislocationAnalysis();

    void setInputCrystalStructure(LatticeStructureType structure);
    void setStructureIdentificationOnly(bool structureIdentificationOnly);
    
    void setMaxTrialCircuitSize(double size);
    void setCircuitStretchability(double stretch);
    void setLineSmoothingLevel(double lineSmoothingLevel);
    void setLinePointInterval(double linePointInterval);
    void setDefectMeshSmoothingLevel(double defectMeshSmoothingLevel);
    
    void setOnlyPerfectDislocations(bool flag);
    void setMarkCoreAtoms(bool markCoreAtoms);
    
    void setIdentificationMode(StructureAnalysis::Mode identificationMode);
    void setRmsd(float rmsd);
    
    json compute(const LammpsParser::Frame &frame, const std::string& jsonOutputFile = "");

private:
    LatticeStructureType _inputCrystalStructure;

    double _maxTrialCircuitSize;
    double _circuitStretchability;
    double _lineSmoothingLevel;
    double _linePointInterval;
    double _defectMeshSmoothingLevel;

    float _rmsd;

    StructureAnalysis::Mode _identificationMode;

    bool _markCoreAtoms;
    bool _structureIdentificationOnly;
    bool _onlyPerfectDislocations;
    
    mutable DXAJsonExporter _jsonExporter;
};

}
