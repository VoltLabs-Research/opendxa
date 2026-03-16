#include <volt/cli/common.h>
#include <volt/core/dislocation_analysis.h>
#include <volt/structures/crystal_structure_types.h>
#include <oneapi/tbb/info.h>
#include <algorithm>
#include <fstream>
#include <set>
#include <string>

using namespace Volt;
using namespace Volt::CLI;

// Local helper functions for parsing
LatticeStructureType parseCrystalStructure(const std::string& str) {
    if (str == "FCC") return LATTICE_FCC;
    if (str == "BCC") return LATTICE_BCC;
    if (str == "HCP") return LATTICE_HCP;
    if (str == "SC")  return LATTICE_SC;
    if (str == "CUBIC_DIAMOND") return LATTICE_CUBIC_DIAMOND;
    if (str == "HEX_DIAMOND")   return LATTICE_HEX_DIAMOND;
    spdlog::warn("Unknown crystal structure '{}', defaulting to BCC.", str);
    return LATTICE_BCC;
}

StructureAnalysis::Mode parseIdentificationMode(const std::string& str) {
    if (str == "CNA") return StructureAnalysis::Mode::CNA;
    if (str == "PTM") return StructureAnalysis::Mode::PTM;
    if (str == "DIAMOND") return StructureAnalysis::Mode::DIAMOND;
    spdlog::warn("Unknown identification mode '{}', defaulting to CNA.", str);
    return StructureAnalysis::Mode::CNA;
}

void showUsage(const std::string& name) {
    printUsageHeader(name, "Volt - Full Dislocation Analysis");
    std::cerr
        << "  --crystalStructure <type>         Reference crystal structure. (BCC|FCC|HCP|CUBIC_DIAMOND|HEX_DIAMOND|SC) [default: BCC]\n"
        << "  --identificationMode <mode>       Structure identification mode. (CNA|PTM|DIAMOND) [default: CNA]\n"
        << "  --rmsd <float>                    RMSD threshold for PTM. [default: 0.1]\n"
        << "  --maxTrialCircuitSize <int>       Maximum Burgers circuit size. [default: 14]\n"
        << "  --circuitStretchability <int>     Circuit stretchability factor. [default: 9]\n"
        << "  --lineSmoothingLevel <float>      Line smoothing level. [default: 1]\n"
        << "  --linePointInterval <float>       Point interval on dislocation lines. [default: 2.5]\n"
        << "  --onlyPerfectDislocations <bool>  Detect only perfect dislocations. [default: false]\n"
        << "  --markCoreAtoms <bool>            Mark dislocation core atoms. [default: false]\n"
        << "  --threads <int>                   Max worker threads (TBB/OMP). [default: auto capped to physical cores]\n";
    printHelpOption();
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        showUsage(argv[0]);
        return 1;
    }
    
    std::string filename, outputBase;
    auto opts = parseArgs(argc, argv, filename, outputBase);
    
    if (hasOption(opts, "--help") || filename.empty()) {
        showUsage(argv[0]);
        return filename.empty() ? 1 : 0;
    }

    if(!hasOption(opts, "--threads")){
        const int maxAvailableThreads = static_cast<int>(oneapi::tbb::info::default_concurrency());
        int physicalCores = 0;
        std::ifstream cpuinfo("/proc/cpuinfo");
        if(cpuinfo.is_open()){
            std::set<std::pair<int, int>> physicalCoreIds;
            int fallbackCpuCores = 0;
            int physicalId = -1;
            int coreId = -1;
            std::string line;

            while(std::getline(cpuinfo, line)){
                if(line.empty()){
                    if(physicalId >= 0 && coreId >= 0){
                        physicalCoreIds.emplace(physicalId, coreId);
                    }
                    physicalId = -1;
                    coreId = -1;
                    continue;
                }

                int parsedValue = 0;
                const auto separator = line.find(':');
                if(separator != std::string::npos){
                    try{
                        parsedValue = std::stoi(line.substr(separator + 1));
                    }catch(const std::exception&){
                        parsedValue = 0;
                    }
                }

                if(line.rfind("physical id", 0) == 0){
                    physicalId = parsedValue;
                }else if(line.rfind("core id", 0) == 0){
                    coreId = parsedValue;
                }else if(line.rfind("cpu cores", 0) == 0){
                    fallbackCpuCores = std::max(fallbackCpuCores, parsedValue);
                }
            }

            if(physicalId >= 0 && coreId >= 0){
                physicalCoreIds.emplace(physicalId, coreId);
            }

            if(!physicalCoreIds.empty()){
                physicalCores = static_cast<int>(physicalCoreIds.size());
            }else{
                physicalCores = fallbackCpuCores;
            }
        }

        const int defaultThreads = std::max(
            1,
            std::min(maxAvailableThreads, physicalCores > 0 ? physicalCores : maxAvailableThreads)
        );
        opts["--threads"] = std::to_string(defaultThreads);
    }
    
    auto parallel = initParallelism(opts, false);
    initLogging("volt-dxa", parallel.threads);
    
    LammpsParser::Frame frame;
    if (!parseFrame(filename, frame)) return 1;
    
    outputBase = deriveOutputBase(filename, outputBase);
    spdlog::info("Output base: {}", outputBase);
    
    DislocationAnalysis analyzer;
    
    analyzer.setInputCrystalStructure(parseCrystalStructure(getString(opts, "--crystalStructure", "BCC")));
    analyzer.setIdentificationMode(parseIdentificationMode(getString(opts, "--identificationMode", "CNA")));
    analyzer.setRmsd(getDouble(opts, "--rmsd", 0.1f));
    analyzer.setMaxTrialCircuitSize(getInt(opts, "--maxTrialCircuitSize", 14));
    analyzer.setCircuitStretchability(getInt(opts, "--circuitStretchability", 9));
    analyzer.setLineSmoothingLevel(getDouble(opts, "--lineSmoothingLevel", 1.0));
    analyzer.setLinePointInterval(getDouble(opts, "--linePointInterval", 2.5));
    analyzer.setOnlyPerfectDislocations(getBool(opts, "--onlyPerfectDislocations"));
    analyzer.setMarkCoreAtoms(getBool(opts, "--markCoreAtoms"));
    
    spdlog::info("Starting dislocation analysis...");
    json result = analyzer.compute(frame, outputBase);
    
    if (result.value("is_failed", false)) {
        spdlog::error("Analysis failed: {}", result.value("error", "Unknown error"));
        return 1;
    }
    
    spdlog::info("Analysis completed successfully.");
    return 0;
}
