#include <volt/cli/common.h>
#include <volt/pipeline/dislocation_analysis.h>
#include <volt/structures/crystal_topology_registry.h>

#include <spdlog/spdlog.h>

#include <exception>
#include <string>

using namespace Volt;
using namespace Volt::CLI;

void showUsage(const std::string& name) {
    printUsageHeader(name, "Volt - Full Dislocation Analysis");
    std::cerr
        << "  --clusters-table <path>           Path to *_clusters.table exported by CNA/PTM.\n"
        << "  --clusters-transitions <path>     Path to *_cluster_transitions.table exported by CNA/PTM.\n"
        << "  --reference-topology <name>       Topology name/alias from OpenDXA YAML definitions for the matrix phase.\n"
        << "  --lattice-dir <path>              Directory containing OpenDXA lattice YAMLs.\n"
        << "  --maxTrialCircuitSize <int>       Maximum Burgers circuit size. [default: 14]\n"
        << "  --circuitStretchability <int>     Circuit stretchability factor. [default: 9]\n"
        << "  --lineSmoothingLevel <float>      Line smoothing level. [default: 1]\n"
        << "  --linePointInterval <float>       Point interval on dislocation lines. [default: 2.5]\n"
        << "  --ghost-layer-scale <float>       Ghost-layer scale relative to max neighbor distance. [default: 3.5]\n"
        << "  --interface-alpha-scale <float>   Interface alpha scale. [default: 5.0]\n"
        << "  --inteface-alpha-scale <float>    Alias for --interface-alpha-scale.\n"
        << "  --crystal-path-steps <int>        Maximum crystal path search depth. [default: 4]\n"
        << "  --export-defect-mesh <bool>       Export defect mesh msgpack. [default: true]\n"
        << "  --export-interface-mesh <bool>    Export interface mesh msgpack. [default: false]\n"
        << "  --export-dislocations <bool>      Export dislocations msgpack. [default: true]\n"
        << "  --export-circuit-information <bool> Export circuit statistics in dislocations msgpack. [default: true]\n"
        << "  --export-dislocation-network-stats <bool> Export network statistics in dislocations msgpack. [default: true]\n"
        << "  --export-junctions <bool>         Export junction information in dislocations msgpack. [default: true]\n"
        << "  --clip-pbc-segments <bool>        Clip dislocation segments against PBC when exporting. [default: true]\n"
        << "  --cover-domain-with-finite-tets <bool> Cover the domain with helper finite tetrahedra. [default: false]\n";
    printHelpOption();
}

int main(int argc, char* argv[]) {
    std::string filename, outputBase;
    auto opts = parseArgs(argc, argv, filename, outputBase);

    if(const int startupStatus = handleHelpOrMissingInput(argc, argv, opts, filename, showUsage);
       startupStatus >= 0){
        return startupStatus;
    }

    initLogging("volt-dxa");
        
    LammpsParser::Frame frame;
    if (!parseFrame(filename, frame)) return 1;
    
    outputBase = deriveOutputBase(filename, outputBase);
    spdlog::info("Output base: {}", outputBase);
    
    DislocationAnalysis analyzer;
    const std::string latticeDirectory = getString(opts, "--lattice-dir", "");
    if(!latticeDirectory.empty()){
        setCrystalTopologySearchRoot(latticeDirectory);
        spdlog::info("Using lattice directory: {}", latticeDirectory);
    }
    if(!hasOption(opts, "--reference-topology")){
        spdlog::error("Missing required option --reference-topology");
        return 1;
    }
    const std::string topologyName = getString(opts, "--reference-topology");
    const CrystalTopologyEntry* topology = crystalTopologyByName(topologyName);
    if(!topology){
        spdlog::error("Unknown value for --reference-topology: {}", topologyName);
        return 1;
    }

    analyzer.setReferenceTopology(topology->name);
    analyzer.setMaxTrialCircuitSize(getInt(opts, "--maxTrialCircuitSize", 14));
    analyzer.setCircuitStretchability(getInt(opts, "--circuitStretchability", 9));
    analyzer.setLineSmoothingLevel(getDouble(opts, "--lineSmoothingLevel", 1.0));
    analyzer.setLinePointInterval(getDouble(opts, "--linePointInterval", 2.5));
    analyzer.setGhostLayerScale(getDouble(opts, "--ghost-layer-scale", 3.5));
    analyzer.setInterfaceAlphaScale(getDouble(
        opts,
        "--interface-alpha-scale",
        getDouble(opts, "--inteface-alpha-scale", 5.0)
    ));
    analyzer.setCrystalPathSteps(getInt(opts, "--crystal-path-steps", 4));
    analyzer.setExportDefectMesh(getBool(opts, "--export-defect-mesh", true));
    analyzer.setExportInterfaceMesh(getBool(opts, "--export-interface-mesh", false));
    analyzer.setExportDislocations(getBool(opts, "--export-dislocations", true));
    analyzer.setExportCircuitInformation(getBool(opts, "--export-circuit-information", true));
    analyzer.setExportDislocationNetworkStats(getBool(opts, "--export-dislocation-network-stats", true));
    analyzer.setExportJunctions(getBool(opts, "--export-junctions", true));
    analyzer.setClipPbcSegments(getBool(opts, "--clip-pbc-segments", true));
    analyzer.setCoverDomainWithFiniteTets(getBool(opts, "--cover-domain-with-finite-tets", false));
    analyzer.setClustersTablePath(getString(opts, "--clusters-table"));
    analyzer.setClusterTransitionsPath(
        getString(opts, "--clusters-transitions", getString(opts, "--cluster-transitions"))
    );
    
    spdlog::info("Starting dislocation analysis...");
    try{
        analyzer.compute(frame, outputBase);
    }catch(const std::exception& error){
        spdlog::error("Analysis failed: {}", error.what());
        return 1;
    }
    
    spdlog::info("Analysis completed successfully.");
    return 0;
}
