#pragma once

#include <volt/core/volt.h>
#include <volt/core/lammps_parser.h>

#include <string>
#include <string_view>
#include <utility>

namespace Volt {

class DislocationAnalysis{
public:
    DislocationAnalysis();

    void compute(const LammpsParser::Frame &frame, const std::string& jsonOutputFile = "");

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
        _exportDefectMesh = exportDefectMesh;
    }

    void setExportInterfaceMesh(bool exportInterfaceMesh){
        _exportInterfaceMesh = exportInterfaceMesh;
    }

    void setExportDislocations(bool exportDislocations){
        _exportDislocations = exportDislocations;
    }

    void setExportCircuitInformation(bool exportCircuitInformation){
        _exportCircuitInformation = exportCircuitInformation;
    }

    void setExportDislocationNetworkStats(bool exportDislocationNetworkStats){
        _exportDislocationNetworkStats = exportDislocationNetworkStats;
    }

    void setExportJunctions(bool exportJunctions){
        _exportJunctions = exportJunctions;
    }

    void setClipPbcSegments(bool clipPbcSegments){
        _clipPbcSegments = clipPbcSegments;
    }

    void setCoverDomainWithFiniteTets(bool coverDomainWithFiniteTets){
        _coverDomainWithFiniteTets = coverDomainWithFiniteTets;
    }

private:
    std::string _referenceTopologyName;

    double _maxTrialCircuitSize;
    double _circuitStretchability;
    double _lineSmoothingLevel;
    double _linePointInterval;
    double _ghostLayerScale;
    double _interfaceAlphaScale;
    int _crystalPathSteps;

    bool _exportDefectMesh;
    bool _exportInterfaceMesh;
    bool _exportDislocations;
    bool _exportCircuitInformation;
    bool _exportDislocationNetworkStats;
    bool _exportJunctions;
    bool _clipPbcSegments;
    bool _coverDomainWithFiniteTets;

    std::string _clustersTablePath;
    std::string _clusterTransitionsPath;
};

}
