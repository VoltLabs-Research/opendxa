#pragma once

#include <volt/core/lammps_parser.h>
#include <volt/helpers/dxa_serialization.h>
#include <volt/plugin/output_files.h>

#include <array>
#include <string>
#include <string_view>
#include <vector>

namespace Volt{

// The files an analysis writes. Core atoms and the defect mesh are only computed when they are written.
struct DislocationExports{
    bool defectMesh = true;
    bool dislocations = true;
    bool interfaceMesh = false;
    bool tessellation = false;
    bool atoms = false;
    bool coreAtoms = false;
    bool coherentRegions = false;
    DxaSerialization::DislocationsExportOptions dislocationTables;
};

// What the files are written from: where they go, which ones, and what the analysis built.
struct DislocationOutputs{
    const std::string& base;
    const DislocationExports& exports;
    const LammpsParser::Frame& frame;
    const StructureAnalysis& structure;
    const DelaunayTessellation& tessellation;
    const InterfaceMesh& interfaceMesh;
    InterfaceMeshTopology& defectMesh;
    const DislocationNetwork& network;
    const std::vector<int>& coreAtomDislocationIds;

    // What writing the files these exports ask for weighs next to the analysis's other steps.
    static double weight(const DislocationExports& exports);

    // Writes each file the exports ask for, one progress step per file.
    void write() const;

    std::string path(std::string_view name) const{ return outputPath(base, name); }

private:
    static void writeCoreAtoms(
        const std::string& filePath,
        const LammpsParser::Frame& frame,
        const StructureAnalysis& structure,
        const std::vector<int>& dislocationIds
    );

    static const std::array<OutputFile<DislocationExports, DislocationOutputs>, 7> files;
};

}
