#include <volt/pipeline/dislocation_outputs.h>

#include <volt/analysis/structure_identification_export.h>
#include <volt/plugin/output_files.h>
#include <volt/structures/crystal_structure_types.h>
#include <volt/utilities/parquet_atom_writer.h>

namespace Volt{

using namespace Volt::Particles;

void DislocationOutputs::writeCoreAtoms(
    const std::string& filePath,
    const LammpsParser::Frame& frame,
    const StructureAnalysis& structure,
    const std::vector<int>& dislocationIds
){
    std::vector<std::size_t> coreAtoms;
    for(std::size_t i = 0; i < dislocationIds.size(); ++i){
        if(dislocationIds[i] >= 0) coreAtoms.push_back(i);
    }

    LammpsParser::Frame coreFrame;
    coreFrame.timestep = frame.timestep;
    coreFrame.natoms = static_cast<int>(coreAtoms.size());
    coreFrame.simulationCell = frame.simulationCell;
    coreFrame.positions.reserve(coreAtoms.size());
    coreFrame.ids.reserve(coreAtoms.size());
    for(const std::size_t source : coreAtoms){
        coreFrame.positions.push_back(source < frame.positions.size()
            ? frame.positions[source] : Point3::Origin());
        coreFrame.ids.push_back(source < frame.ids.size()
            ? frame.ids[source] : static_cast<int>(source));
    }

    const StructureContext& context = structure.context();
    const int* structureTypes = context.structureTypes ? context.structureTypes->constDataInt() : nullptr;

    streamAtomsToParquet(
        filePath,
        coreFrame,
        [&coreAtoms, structureTypes](std::size_t atomIndex){
            return std::string(structureTypeName(
                structureTypes
                    ? structureTypes[coreAtoms[atomIndex]]
                    : static_cast<int>(StructureType::OTHER)
            ));
        },
        [&coreAtoms, &dislocationIds](ColumnarAtomWriter& writer, std::size_t atomIndex){
            writer.field("dislocation_id", dislocationIds[coreAtoms[atomIndex]]);
        }
    );
}

const std::array<OutputFile<DislocationExports, DislocationOutputs>, 7> DislocationOutputs::files{{
    {"Writing the defect mesh", 1.0, &DislocationExports::defectMesh, [](const DislocationOutputs& out){
        DxaSerialization::streamMeshToFile(out.path("defect_mesh"), out.defectMesh, out.frame.simulationCell);
        out.defectMesh.clear();
    }},
    {"Writing dislocations", 0.2, &DislocationExports::dislocations, [](const DislocationOutputs& out){
        DxaSerialization::streamDislocationsToFile(
            out.path("dislocations"), out.path("dislocation_summary"),
            out.network, out.structure, out.exports.dislocationTables
        );
    }},
    {"Writing the interface mesh", 1.0, &DislocationExports::interfaceMesh, [](const DislocationOutputs& out){
        DxaSerialization::streamInterfaceMeshToFile(out.path("interface_mesh"), out.interfaceMesh);
    }},
    {"Writing the Delaunay tessellation", 8.0, &DislocationExports::tessellation, [](const DislocationOutputs& out){
        DxaSerialization::streamDelaunayTessellationToFile(out.path("delaunay_tessellation"), out.tessellation);
    }},
    {"Writing the structure identification", 2.0, &DislocationExports::atoms, [](const DislocationOutputs& out){
        StructureIdentificationExport::streamStructureIdentificationToParquet(out.path("atoms"), out.frame, out.structure);
    }},
    {"Writing core atoms", 0.2, &DislocationExports::coreAtoms, [](const DislocationOutputs& out){
        writeCoreAtoms(out.path("core_atoms"), out.frame, out.structure, out.coreAtomDislocationIds);
    }},
    {"Writing coherent crystalline regions", 3.0, &DislocationExports::coherentRegions, [](const DislocationOutputs& out){
        DxaSerialization::streamCoherentCrystallineRegionsToFile(
            out.path("coherent_crystalline_regions"), out.frame, out.structure
        );
    }},
}};

double DislocationOutputs::weight(const DislocationExports& exports){
    return totalWeight(outputSteps(files, exports));
}

void DislocationOutputs::write() const{
    writeOutputFiles(files, exports, *this);
}

}
