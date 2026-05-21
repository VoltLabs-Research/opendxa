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
        << "  --clusters_table <path>           Path to *_clusters.table exported by CNA/PTM.\n"
        << "  --clusters_transitions <path>     Path to *_cluster_transitions.table exported by CNA/PTM.\n"
        << "  --reference_topology <name>       Topology name/alias from OpenDXA YAML definitions for the matrix phase.\n"
        << "  --lattice_dir <path>              Directory containing OpenDXA lattice YAMLs.\n"
        << "  --max_trial_circuit_size <int>    Maximum Burgers circuit size. [default: 14]\n"
        << "  --circuit_stretchability <int>    Circuit stretchability factor. [default: 9]\n"
        << "  --line_smoothing_level <float>    Line smoothing level. [default: 1]\n"
        << "  --line_point_interval <float>     Point interval on dislocation lines. [default: 2.5]\n"
        << "  --ghost_layer_scale <float>       Ghost-layer scale relative to max neighbor distance. [default: 3.5]\n"
        << "  --interface_alpha_scale <float>   Interface alpha scale. [default: 5.0]\n"
        << "  --inteface_alpha_scale <float>    Accepted alias for --interface_alpha_scale.\n"
        << "  --crystal_path_steps <int>        Maximum crystal path search depth. [default: 4]\n"
        << "  --export_defect_mesh <bool>       Export defect mesh msgpack. [default: true]\n"
        << "  --export_interface_mesh <bool>    Export interface mesh msgpack. [default: false]\n"
        << "  --export_delaunay_tessellation <bool> Export Delaunay tessellation msgpack. [default: false]\n"
        << "  --export_structure_identification <bool> Export atoms.msgpack grouped by structure id/name. [default: false]\n"
        << "  --export_coherent_crystalline_regions <bool> Export atoms grouped by coherent cluster id. [default: false]\n"
        << "  --export_dislocations <bool>      Export dislocations msgpack. [default: true]\n"
        << "  --export_circuit_information <bool> Export circuit statistics in dislocations msgpack. [default: true]\n"
        << "  --export_dislocation_network_stats <bool> Export network statistics in dislocations msgpack. [default: true]\n"
        << "  --export_junctions <bool>         Export junction information in dislocations msgpack. [default: true]\n"
        << "  --clip_pbc_segments <bool>        Clip dislocation segments against PBC when exporting. [default: true]\n"
        << "  --cover_domain_with_finite_tets <bool> Cover the domain with helper finite tetrahedra. [default: false]\n";
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
    const std::string latticeDirectory = getString(opts, "--lattice_dir", "");
    if(!latticeDirectory.empty()){
        setCrystalTopologySearchRoot(latticeDirectory);
        spdlog::info("Using lattice directory: {}", latticeDirectory);
    }
    if(!hasOption(opts, "--reference_topology")){
        spdlog::error("Missing required option --reference_topology");
        return 1;
    }
    const std::string topologyName = getString(opts, "--reference_topology");
    const CrystalTopologyEntry* topology = crystalTopologyByName(topologyName);
    if(!topology){
        spdlog::error("Unknown value for --reference_topology: {}", topologyName);
        return 1;
    }

    analyzer.setReferenceTopology(topology->name);
    analyzer.setMaxTrialCircuitSize(getInt(opts, "--max_trial_circuit_size", 14));
    analyzer.setCircuitStretchability(getInt(opts, "--circuit_stretchability", 9));
    analyzer.setLineSmoothingLevel(getDouble(opts, "--line_smoothing_level", 1.0));
    analyzer.setLinePointInterval(getDouble(opts, "--line_point_interval", 2.5));
    analyzer.setGhostLayerScale(getDouble(opts, "--ghost_layer_scale", 3.5));
    analyzer.setInterfaceAlphaScale(getDouble(
        opts,
        "--interface_alpha_scale",
        getDouble(opts, "--inteface_alpha_scale", 5.0)
    ));
    analyzer.setCrystalPathSteps(getInt(opts, "--crystal_path_steps", 4));
    analyzer.setExportDefectMesh(getBool(opts, "--export_defect_mesh", true));
    analyzer.setExportInterfaceMesh(getBool(opts, "--export_interface_mesh", false));
    analyzer.setExportDelaunayTessellation(getBool(opts, "--export_delaunay_tessellation", false));
    analyzer.setExportStructureIdentification(getBool(opts, "--export_structure_identification", false));
    analyzer.setExportCoherentCrystallineRegions(getBool(opts, "--export_coherent_crystalline_regions", false));
    analyzer.setExportDislocations(getBool(opts, "--export_dislocations", true));
    analyzer.setExportCircuitInformation(getBool(opts, "--export_circuit_information", true));
    analyzer.setExportDislocationNetworkStats(getBool(opts, "--export_dislocation_network_stats", true));
    analyzer.setExportJunctions(getBool(opts, "--export_junctions", true));
    analyzer.setClipPbcSegments(getBool(opts, "--clip_pbc_segments", true));
    analyzer.setCoverDomainWithFiniteTets(getBool(opts, "--cover_domain_with_finite_tets", false));
    analyzer.setClustersTablePath(getString(opts, "--clusters_table"));
    analyzer.setClusterTransitionsPath(getString(opts, "--clusters_transitions"));
    
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
