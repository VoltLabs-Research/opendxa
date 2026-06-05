#include <volt/plugin/plugin_entry.h>
#include <volt/pipeline/dislocation_analysis.h>
#include <volt/structures/crystal_topology_registry.h>

using namespace Volt;
using namespace Volt::Plugin;
using S = DislocationAnalysis;

static const std::vector<OptionBinding<S>> bindings = {
    opt("--clusters_table", "Clusters table path", "", &S::setClustersTablePath),
    opt("--clusters_transitions", "Transitions path", "", &S::setClusterTransitionsPath),
    opt("--max_trial_circuit_size", "Max trial circuit size", 14.0, &S::setMaxTrialCircuitSize),
    opt("--circuit_stretchability", "Circuit stretchability", 9.0, &S::setCircuitStretchability),
    opt("--line_smoothing_level", "Line smoothing level", 1.0, &S::setLineSmoothingLevel),
    opt("--line_point_interval", "Point interval on dislocation lines", 2.5, &S::setLinePointInterval),
    opt("--ghost_layer_scale", "Ghost layer scale", 3.5, &S::setGhostLayerScale),
    opt("--interface_alpha_scale", "Interface alpha scale", 5.0, &S::setInterfaceAlphaScale),
    opt("--crystal_path_steps", "Crystal path search depth", 4, &S::setCrystalPathSteps),
    opt("--export_defect_mesh", "Export defect mesh", true, &S::setExportDefectMesh),
    opt("--export_interface_mesh", "Export interface mesh", false, &S::setExportInterfaceMesh),
    opt("--export_delaunay_tessellation", "Export Delaunay tessellation", false, &S::setExportDelaunayTessellation),
    opt("--export_structure_identification", "Export structure identification", false, &S::setExportStructureIdentification),
    opt("--export_coherent_crystalline_regions", "Export coherent crystalline regions", false, &S::setExportCoherentCrystallineRegions),
    opt("--export_dislocations", "Export dislocations", true, &S::setExportDislocations),
    opt("--export_circuit_information", "Export circuit information", true, &S::setExportCircuitInformation),
    opt("--export_dislocation_network_stats", "Export network stats", true, &S::setExportDislocationNetworkStats),
    opt("--export_junctions", "Export junctions", true, &S::setExportJunctions),
    opt("--clip_pbc_segments", "Clip PBC segments", true, &S::setClipPbcSegments),
    opt("--cover_domain_with_finite_tets", "Cover domain with finite tets", false, &S::setCoverDomainWithFiniteTets),
};

VOLT_PLUGIN_MAIN(
    (PluginDescriptor{"volt-dxa", "Full Dislocation Analysis", optionsMeta(bindings)}),
    [](const OptsMap& opts, const LammpsParser::Frame& frame,
       const LammpsParser::Frame*, const std::string& outputBase) -> json {
    const std::string latticeDir = CLI::getString(opts, "--lattice_dir", "");
    if (!latticeDir.empty()) setCrystalTopologySearchRoot(latticeDir);

    if (auto err = requireOptions(opts, {"--reference_topology"})) return *err;

    const std::string topoName = CLI::getString(opts, "--reference_topology");
    const auto* topo = crystalTopologyByName(topoName);
    if (!topo)
        return AnalysisResult::failure("Unknown reference topology: " + topoName);

    S analyzer;
    analyzer.setReferenceTopology(topo->name);
    applyAll(analyzer, bindings, opts);

    try { analyzer.compute(frame, outputBase); }
    catch (const std::exception& e) {
        return AnalysisResult::failure(std::string("DXA failed: ") + e.what());
    }

    json result;
    result["is_failed"] = false;
    return result;
})
