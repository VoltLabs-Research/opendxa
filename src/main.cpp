#include <volt/cli/common.h>
#include <volt/pipeline/dislocation_analysis.h>

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
        << "  --reference-structure-label <int> Structure label in *_clusters.table that represents the reference matrix phase.\n"
        << "  --maxTrialCircuitSize <int>       Maximum Burgers circuit size. [default: 14]\n"
        << "  --circuitStretchability <int>     Circuit stretchability factor. [default: 9]\n"
        << "  --lineSmoothingLevel <float>      Line smoothing level. [default: 1]\n"
        << "  --linePointInterval <float>       Point interval on dislocation lines. [default: 2.5]\n";
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
    if(!hasOption(opts, "--reference-structure-label")){
        spdlog::error("Missing required option --reference-structure-label");
        return 1;
    }
    int referenceStructureLabel = 0;
    try{
        referenceStructureLabel = std::stoi(getString(opts, "--reference-structure-label"));
    }catch(...){
        spdlog::error("Invalid value for --reference-structure-label");
        return 1;
    }

    analyzer.setReferenceStructureLabel(referenceStructureLabel);
    analyzer.setMaxTrialCircuitSize(getInt(opts, "--maxTrialCircuitSize", 14));
    analyzer.setCircuitStretchability(getInt(opts, "--circuitStretchability", 9));
    analyzer.setLineSmoothingLevel(getDouble(opts, "--lineSmoothingLevel", 1.0));
    analyzer.setLinePointInterval(getDouble(opts, "--linePointInterval", 2.5));
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
