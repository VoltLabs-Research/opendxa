#pragma once

#include <volt/core/volt.h>
#include <volt/core/lammps_parser.h>
#include <volt/pipeline/dislocation_outputs.h>

#include <string>
#include <utility>

namespace Volt {

struct ClusterGraphExportPaths;
class ReconstructedStructureContext;

class DislocationAnalysis{
public:
    void compute(const LammpsParser::Frame& frame, const std::string& outputBase);

    void setReferenceTopology(std::string topologyName){
        _referenceTopologyName = std::move(topologyName);
    }

    void setMaxTrialCircuitSize(double size){
        _maxTrialCircuitSize = size;
    }

    void setCircuitStretchability(double stretch){
        _circuitStretchability = stretch;
    }

    void setClustersTablePath(std::string path){
        _clustersTablePath = std::move(path);
    }

    void setClusterTransitionsPath(std::string path){
        _clusterTransitionsPath = std::move(path);
    }

    void setNeighborLatticePath(std::string path){
        _neighborLatticePath = std::move(path);
    }

    void setLineSmoothingLevel(double lineSmoothingLevel){
        _lineSmoothingLevel = lineSmoothingLevel;
    }

    void setLinePointInterval(double linePointInterval){
        _linePointInterval = linePointInterval;
    }

    void setGhostLayerScale(double ghostLayerScale){
        _ghostLayerScale = ghostLayerScale;
    }

    void setInterfaceAlphaScale(double interfaceAlphaScale){
        _interfaceAlphaScale = interfaceAlphaScale;
    }

    void setCrystalPathSteps(int crystalPathSteps){
        _crystalPathSteps = crystalPathSteps;
    }

    void setExportDefectMesh(bool exportDefectMesh){
        _exports.defectMesh = exportDefectMesh;
    }

    void setExportInterfaceMesh(bool exportInterfaceMesh){
        _exports.interfaceMesh = exportInterfaceMesh;
    }

    void setExportDelaunayTessellation(bool exportDelaunayTessellation){
        _exports.tessellation = exportDelaunayTessellation;
    }

    void setExportStructureIdentification(bool exportStructureIdentification){
        _exports.atoms = exportStructureIdentification;
    }

    void setMarkCoreAtoms(bool markCoreAtoms){
        _exports.coreAtoms = markCoreAtoms;
    }

    void setExportCoherentCrystallineRegions(bool exportCoherentCrystallineRegions){
        _exports.coherentRegions = exportCoherentCrystallineRegions;
    }

    void setExportDislocations(bool exportDislocations){
        _exports.dislocations = exportDislocations;
    }

    void setExportCircuitInformation(bool exportCircuitInformation){
        _exports.dislocationTables.exportCircuitInformation = exportCircuitInformation;
    }

    void setExportDislocationNetworkStats(bool exportDislocationNetworkStats){
        _exports.dislocationTables.exportDislocationNetworkStats = exportDislocationNetworkStats;
    }

    void setExportJunctions(bool exportJunctions){
        _exports.dislocationTables.exportJunctions = exportJunctions;
    }

    void setClipPbcSegments(bool clipPbcSegments){
        _exports.dislocationTables.clipPbcSegments = clipPbcSegments;
    }

    void setCoverDomainWithFiniteTets(bool coverDomainWithFiniteTets){
        _coverDomainWithFiniteTets = coverDomainWithFiniteTets;
    }

private:
    static void reconstructStructure(
        const LammpsParser::Frame& frame,
        const std::string& neighborLatticePath,
        const ClusterGraphExportPaths& clusterPaths,
        StructureAnalysis& structure,
        ReconstructedStructureContext& context
    );
    static void mapElastically(ElasticMapping& elasticMap, int crystalPathSteps);

    std::string _referenceTopologyName;
    std::string _clustersTablePath;
    std::string _clusterTransitionsPath;
    std::string _neighborLatticePath;

    double _maxTrialCircuitSize = 14.0;
    double _circuitStretchability = 9.0;
    double _lineSmoothingLevel = 1.0;
    double _linePointInterval = 2.5;
    double _ghostLayerScale = 3.5;
    double _interfaceAlphaScale = 5.0;
    int _crystalPathSteps = 4;
    bool _coverDomainWithFiniteTets = false;

    DislocationExports _exports;
};

}
