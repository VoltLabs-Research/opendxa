#pragma once

#include <volt/core/lammps_parser.h>
#include <volt/core/run_log.h>
#include <volt/pipeline/interface_mesh.h>

#include <cstddef>
#include <vector>

namespace Volt{

class StructureAnalysis;
class StructureContext;
class DelaunayTessellation;
class ElasticMapping;
class BurgersLoopBuilder;
class DislocationNetwork;
struct DislocationExports;

// The execution log of one analysis. Each phase hands it what it built, and the log writes what that holds and how
// long the phase took. Scripts parse these lines, so their wording stays as it is.
class DislocationAnalysisLog : public RunLog{
public:
    explicit DislocationAnalysisLog(const LammpsParser::Frame& frame);

    using RunLog::dump;

    void dump(const StructureAnalysis& structure);
    void dump(const DelaunayTessellation& tessellation);
    void dump(const ElasticMapping& elasticMap);
    void dump(const InterfaceMesh& interfaceMesh);
    void dump(const BurgersLoopBuilder& tracer);
    void dump(const std::vector<int>& coreAtomDislocationIds);
    void dump(const InterfaceMeshTopology& defectMesh);
    void dump(const DislocationNetwork& smoothedNetwork);
    // The export is the last phase, so its line closes the run with the total.
    void dump(const DislocationExports& exports);

private:
    struct NetworkTopology{
        std::size_t closedLoops = 0;
        std::size_t infiniteLines = 0;
        std::size_t junctions = 0;
        std::size_t danglingEnds = 0;
    };

    static NetworkTopology summarizeNetworkTopology(const DislocationNetwork& network);
    static std::size_t clusteredAtomCount(const StructureContext& context);

    const double _cellVolume;
    std::size_t _interfaceFacets = 0;
};

}
