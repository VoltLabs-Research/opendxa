#include <volt/plugin/plugin_entry.h>
#include <volt/pipeline/dislocation_analysis.h>
#include <volt/structures/crystal_topology_registry.h>

#include <algorithm>
#include <array>
#include <string>
#include <vector>

using namespace Volt;
using namespace Volt::Plugin;
using S = DislocationAnalysis;

static const std::vector<OptionBinding<S>> bindings = {
    opt("--clusters_table", "Clusters table path", "", &S::setClustersTablePath),
    opt("--clusters_transitions", "Transitions path", "", &S::setClusterTransitionsPath),
    opt("--neighbor_lattice", "Per-atom neighbor topology parquet path", "", &S::setNeighborLatticePath),
    opt("--max_trial_circuit_size", "Max trial circuit size", 14.0, &S::setMaxTrialCircuitSize),
    opt("--circuit_stretchability", "Circuit stretchability", 9.0, &S::setCircuitStretchability),
    opt("--line_smoothing_level", "Line smoothing level", 1.0, &S::setLineSmoothingLevel),
    opt("--line_point_interval", "Point interval on dislocation lines", 2.5, &S::setLinePointInterval),
    opt("--ghost_layer_scale", "Ghost layer scale", 3.5, &S::setGhostLayerScale),
    opt("--interface_alpha_scale", "Interface alpha scale", 5.0, &S::setInterfaceAlphaScale),
    opt("--crystal_path_steps", "Crystal path search depth", 4, &S::setCrystalPathSteps),
    opt("--export_defect_mesh", "Export defect mesh", true, &S::setExportDefectMesh),
    opt("--export_interface_mesh", "Export interface mesh (currently byte-identical to the defect mesh)", false, &S::setExportInterfaceMesh),
    opt("--export_delaunay_tessellation", "Export Delaunay tessellation", false, &S::setExportDelaunayTessellation),
    opt("--export_structure_identification", "Export per-atom structure identification", false, &S::setExportStructureIdentification),
    opt("--mark_core_atoms", "Mark dislocation core atoms (adds a dislocation_id column to the per-atom table)", false, &S::setMarkCoreAtoms),
    opt("--export_coherent_crystalline_regions", "Export coherent crystalline regions", false, &S::setExportCoherentCrystallineRegions),
    opt("--export_dislocations", "Export dislocations", true, &S::setExportDislocations),   
    opt("--export_circuit_information", "Export circuit information", true, &S::setExportCircuitInformation),
    opt("--export_dislocation_network_stats", "Export network stats", true, &S::setExportDislocationNetworkStats),
    opt("--export_junctions", "Export junctions", true, &S::setExportJunctions),
    opt("--clip_pbc_segments", "Clip PBC segments", true, &S::setClipPbcSegments),
    opt("--cover_domain_with_finite_tets", "Cover domain with finite tets", false, &S::setCoverDomainWithFiniteTets),
};

static void refineTopologyDomain(std::vector<CliOption>& options, const OptsMap& opts) {
    const auto latticeDir = CLI::getString(opts, "--lattice_dir", "");
    if(latticeDir.empty()){
        return;
    }
    std::vector<std::string> names;
    try{
        setCrystalTopologySearchRoot(latticeDir);
        for(const auto& entry : crystalTopologyRegistry().entries()){
            names.push_back(entry.name);
        }
    }catch(const std::exception&){
        return;
    }
    if(names.empty()){
        return;
    }
    std::sort(names.begin(), names.end());
    for(auto& option : options){
        if(option.name == "--reference_topology"){
            option.type = "enum";
            option.values = std::move(names);
            return;
        }
    }
}

static PluginDescriptor buildDescriptor() {
    auto meta = optionsMeta(bindings);
    meta.insert(meta.begin(), {"--lattice_dir", "path",
        "Directory containing reference lattice/topology YAMLs", "", {}, "share/volt/lattices"});
    meta.insert(meta.begin() + 1, {"--reference_topology", "string",
        "Reference structure name the Burgers vectors are expressed in. OpenDXA consumes the "
        "contract (--clusters_table + --clusters_transitions + --neighbor_lattice) produced by an "
        "upstream structure-identification step (e.g. PolyhedralTemplateMatching, "
        "PatternStructureMatching or ReferenceCrystalMatching). If the name matches a known topology "
        "YAML in --lattice_dir it is canonicalized; otherwise the contract-supplied name is used "
        "as-is. [required]", ""});
    meta.insert(meta.begin() + 2, {"--metric_rescale", "sx,sy,sz",
        "Affine metric isotropization factors for anisotropic crystals. DXA runs in a frame with "
        "coords/(sx,sy,sz); ideal vectors are rescaled to match, and the Burgers vector comes out in "
        "lattice units (multiply components by sx,sy,sz for Angstrom). For reference-built contracts, "
        "use the factors reported by the ReferenceCrystalMatching producer. Default: off (1,1,1).", ""});
    PluginDescriptor descriptor{"volt-dxa", "Full Dislocation Analysis", std::move(meta)};
    descriptor.refine = refineTopologyDomain;
    return descriptor;
}

static bool parseMetricRescale(const std::string& raw, double& rescaleX, double& rescaleY, double& rescaleZ){
    if(raw.empty()){
        return false;
    }
    std::array<double, 3> factors{1.0, 1.0, 1.0};
    std::size_t start = 0;
    for(int axis = 0; axis < 3; ++axis){
        const std::size_t comma = raw.find(',', start);
        const std::string token = raw.substr(start, comma == std::string::npos ? std::string::npos : comma - start);
        try{
            factors[axis] = std::stod(token);
        }catch(...){
            return false;
        }
        if(factors[axis] <= 0.0){
            return false;
        }
        if(comma == std::string::npos){
            if(axis != 2){
                return false;
            }
            break;
        }
        start = comma + 1;
    }
    rescaleX = factors[0];
    rescaleY = factors[1];
    rescaleZ = factors[2];
    return true;
}

VOLT_PLUGIN_MAIN(buildDescriptor(),
    [](const OptsMap& opts, const LammpsParser::Frame& frame,
       const LammpsParser::Frame*, const std::string& outputBase) -> json {
    if(auto err = requireOptions(opts, {"--reference_topology"})){
        return *err;
    }

    const auto latticeDir = CLI::getString(opts, "--lattice_dir", "");
    if(!latticeDir.empty()){
        setCrystalTopologySearchRoot(latticeDir);
    }

    const auto referenceName = CLI::getString(opts, "--reference_topology");

    S analyzer;
    applyAll(analyzer, bindings, opts);

    std::string resolvedTopology = referenceName;
    try{
        if(const auto* topology = crystalTopologyByName(referenceName)){
            resolvedTopology = topology->name;
        }
    }catch(const std::exception&){
    }
    analyzer.setReferenceTopology(resolvedTopology);

    const auto metricRescaleRaw = CLI::getString(opts, "--metric_rescale", "");
    if(!metricRescaleRaw.empty()){
        double rescaleX = 1.0, rescaleY = 1.0, rescaleZ = 1.0;
        if(!parseMetricRescale(metricRescaleRaw, rescaleX, rescaleY, rescaleZ)){
            return AnalysisResult::failure("Malformed --metric_rescale (expected positive sx,sy,sz)");
        }
        analyzer.setMetricRescale(rescaleX, rescaleY, rescaleZ);
    }

    try{
        analyzer.compute(frame, outputBase);
    }catch(const std::exception& e){
        return AnalysisResult::failure(std::string("DXA failed: ") + e.what());
    }

    json result;
    result["is_failed"] = false;
    return result;
})
