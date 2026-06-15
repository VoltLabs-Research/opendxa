#include <volt/plugin/plugin_entry.h>
#include <volt/pipeline/dislocation_analysis.h>
#include <volt/structures/crystal_topology_registry.h>

#include <yaml-cpp/yaml.h>

#include <array>
#include <filesystem>
#include <optional>
#include <string>

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

static PluginDescriptor buildDescriptor() {
    auto meta = optionsMeta(bindings);
    meta.insert(meta.begin(), {"--lattice_dir", "path", "Directory containing lattice/structure YAMLs", ""});
    meta.insert(meta.begin() + 1, {"--reference_topology", "string",
        "Name of the reference structure YAML (in --lattice_dir). A standard topology YAML "
        "(coordination_number + neighbor_vectors) selects classic mode (consumes --clusters_table). "
        "An exotic YAML carrying 'cation_species' + 'reference_crystal' selects full-crystal mode: "
        "OpenDXA builds its own elastic mapping (neighbour graph + ideal vectors + single grain "
        "cluster) from the perfect reference, with auto-selected cutoff and auto metric isotropization. "
        "[required]", ""});
    meta.insert(meta.begin() + 2, {"--metric_rescale", "sx,sy,sz",
        "Classic-mode only: affine metric isotropization factors for anisotropic crystals. "
        "DXA runs in a frame with coords/(sx,sy,sz); Burgers vector comes out in lattice units "
        "(multiply components by sx,sy,sz for Angstrom). Full-crystal mode auto-derives these from "
        "the reference cell. Default: off (1,1,1).", ""});
    return {"volt-dxa", "Full Dislocation Analysis", std::move(meta)};
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

struct ExoticStructure{
    std::string referenceCrystal;
    int cationSpecies = 0;
    double fullCrystalCutoff = 0.0;
};

static std::optional<ExoticStructure> tryLoadExoticStructure(
    const std::string& latticeDir, const std::string& name)
{
    if(latticeDir.empty()){
        return std::nullopt;
    }

    namespace fs = std::filesystem;
    fs::path yamlPath;
    for(const char* extension : {".yml", ".yaml"}){
        const fs::path candidate = fs::path(latticeDir) / (name + extension);
        std::error_code errorCode;
        if(fs::exists(candidate, errorCode)){
            yamlPath = candidate;
            break;
        }
    }
    if(yamlPath.empty()){
        return std::nullopt;
    }

    YAML::Node document;
    try{
        document = YAML::LoadFile(yamlPath.string());
    }catch(const std::exception&){
        return std::nullopt;
    }

    if(!document || !document.IsMap() || !document["cation_species"]){
        return std::nullopt;
    }

    ExoticStructure exotic;
    exotic.cationSpecies = document["cation_species"].as<int>();
    if(!document["reference_crystal"]){
        throw std::runtime_error("Exotic structure '" + name + "' has cation_species but no reference_crystal");
    }

    fs::path referencePath(document["reference_crystal"].as<std::string>());
    if(referencePath.is_relative()){
        referencePath = yamlPath.parent_path() / referencePath;
    }
    std::error_code errorCode;
    if(!fs::exists(referencePath, errorCode)){
        throw std::runtime_error("Exotic structure '" + name + "': reference_crystal not found at '" + referencePath.string() + "'");
    }
    exotic.referenceCrystal = referencePath.string();

    if(document["full_crystal_cutoff"]){
        exotic.fullCrystalCutoff = document["full_crystal_cutoff"].as<double>();
    }
    return exotic;
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

    std::optional<ExoticStructure> exotic;
    try{
        exotic = tryLoadExoticStructure(latticeDir, referenceName);
    }catch(const std::exception& e){
        return AnalysisResult::failure(std::string("Reference structure: ") + e.what());
    }

    if(exotic){
        analyzer.setReferenceTopology(referenceName);
        analyzer.setReferenceCrystal(exotic->referenceCrystal);
        analyzer.setCationSpecies(exotic->cationSpecies);
        analyzer.setFullCrystalCutoff(exotic->fullCrystalCutoff);
    }else{
        const auto* topology = crystalTopologyByName(referenceName);
        if(!topology){
            return AnalysisResult::failure("Unknown --reference_topology '" + referenceName +
                "' (not a known topology, and no exotic YAML with cation_species found in --lattice_dir)");
        }
        analyzer.setReferenceTopology(topology->name);

        const auto metricRescaleRaw = CLI::getString(opts, "--metric_rescale", "");
        if(!metricRescaleRaw.empty()){
            double rescaleX = 1.0, rescaleY = 1.0, rescaleZ = 1.0;
            if(!parseMetricRescale(metricRescaleRaw, rescaleX, rescaleY, rescaleZ)){
                return AnalysisResult::failure("Malformed --metric_rescale (expected positive sx,sy,sz)");
            }
            analyzer.setMetricRescale(rescaleX, rescaleY, rescaleZ);
        }
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
