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
        _exportDefectMesh = exportDefectMesh;
    }

    void setExportInterfaceMesh(bool exportInterfaceMesh){
        _exportInterfaceMesh = exportInterfaceMesh;
    }

    void setExportDelaunayTessellation(bool exportDelaunayTessellation){
        _exportDelaunayTessellation = exportDelaunayTessellation;
    }

    void setExportStructureIdentification(bool exportStructureIdentification){
        _exportStructureIdentification = exportStructureIdentification;
    }

    void setExportCoherentCrystallineRegions(bool exportCoherentCrystallineRegions){
        _exportCoherentCrystallineRegions = exportCoherentCrystallineRegions;
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

    // Affine metric pre-conditioning (isotropization).
    // For strongly ansiotropic crystals (a != b != c), the Delaunay tessellation
    // of the physical point cloud becomes degenerate along the long axis, so no
    // Burgers circuit can close. Setting a per-axis rescale (sa, sb, sc) makes DXA build
    // the tessellation in an isotropized frame (positions / s) where the lattice is
    // (near-)cubic; the ideal lattice vectors that determine the Burgers vector are
    // untouched, so b = sum(deltaX) comes out correct in lattice units. 
    // Exported line geometry is rescale back to physical coordiantes by (sa, sb, sc).
    void setMetricRescale(double sa, double sb, double sc){
        _metricRescaleX = sa;
        _metricRescaleY = sb;
        _metricRescaleZ = sc;
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

    double _metricRescaleX = 1.0;
    double _metricRescaleY = 1.0;
    double _metricRescaleZ = 1.0;

    bool _exportDefectMesh;
    bool _exportInterfaceMesh;
    bool _exportDelaunayTessellation;
    bool _exportStructureIdentification;
    bool _exportCoherentCrystallineRegions;
    bool _exportDislocations;
    bool _exportCircuitInformation;
    bool _exportDislocationNetworkStats;
    bool _exportJunctions;
    bool _clipPbcSegments;
    bool _coverDomainWithFiniteTets;

    std::string _clustersTablePath;
    std::string _clusterTransitionsPath;
    std::string _neighborLatticePath;
};

}
