#include <volt/helpers/full_crystal_context.h>
#include <volt/structures/cluster_graph.h>

#include <ptm_polar.h>
#include <spdlog/spdlog.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <limits>
#include <sstream>
#include <unordered_map>
#include <vector>

namespace Volt{

namespace{

constexpr int kNeighborSlots = MAX_NEIGHBORS;

double medianOf(std::vector<double> values){
    if(values.empty()){
        return 0.0;
    }
    std::sort(values.begin(), values.end());
    return values[values.size() / 2];
}

double percentileOf(std::vector<double> values, double fraction){
    if(values.empty()){
        return 0.0;
    }
    std::sort(values.begin(), values.end());
    const double position = fraction * (static_cast<double>(values.size()) - 1.0);
    const std::size_t lowerIndex = static_cast<std::size_t>(std::floor(position));
    const std::size_t upperIndex = static_cast<std::size_t>(std::ceil(position));
    if(lowerIndex == upperIndex){
        return values[lowerIndex];
    }
    const double weight = position - static_cast<double>(lowerIndex);
    return values[lowerIndex] * (1.0 - weight) + values[upperIndex] * weight;
}

double ambiguityRho(double cellLengthA, double cellLengthB, double cellLengthC){
    return 0.5 * std::min({cellLengthA, cellLengthB, cellLengthC});
}

struct BasisSite{
    int species = 0;
    Vector3 fractional = Vector3::Zero();
    std::vector<std::pair<int, Vector3>> shell;
};

struct PerfectReference{
    bool ok = false;
    std::string message;
    double cellLengthA = 0;
    double cellLengthB = 0;
    double cellLengthC = 0;
    std::vector<BasisSite> sites;
};

struct CationReference{
    double cellLengthA = 0;
    double cellLengthB = 0;
    double cellLengthC = 0;
    std::vector<Vector3> idealFractional;
    int targetSpecies = 0;
    static constexpr double physicalBondCutoff = 5.7;
};

struct GrainFrame{
    bool ok = false;
    std::string message;
    Matrix3 rotation = Matrix3::Identity();
    double bulkSnapResidual = 0;
    int cationAtomCount = 0;
};

struct AllAtomNeighbors{
    std::vector<int> offsets;
    std::vector<int> indices;
    std::vector<Vector3> deltas;
};

PerfectReference buildPerfectReference(const std::string& referenceFile, double cutoff){
    PerfectReference reference;

    std::ifstream input(referenceFile);
    if(!input){
        reference.message = "cannot open reference file '" + referenceFile + "'";
        return reference;
    }

    double boxXLow = 0, boxXHigh = 0;
    double boxYLow = 0, boxYHigh = 0;
    double boxZLow = 0, boxZHigh = 0;
    bool haveX = false, haveY = false, haveZ = false;

    std::vector<int> species;
    std::vector<std::array<double, 3>> positions;
    std::string line;
    bool insideAtomsSection = false;

    while(std::getline(input, line)){
        const std::size_t firstNonSpace = line.find_first_not_of(" \t");
        if(firstNonSpace == std::string::npos){
            continue;
        }
        const std::string trimmed = line.substr(firstNonSpace);

        if(!insideAtomsSection){
            if(trimmed.find("xlo") != std::string::npos){
                std::istringstream stream(trimmed);
                stream >> boxXLow >> boxXHigh;
                haveX = true;
                continue;
            }
            if(trimmed.find("ylo") != std::string::npos){
                std::istringstream stream(trimmed);
                stream >> boxYLow >> boxYHigh;
                haveY = true;
                continue;
            }
            if(trimmed.find("zlo") != std::string::npos){
                std::istringstream stream(trimmed);
                stream >> boxZLow >> boxZHigh;
                haveZ = true;
                continue;
            }
            if(trimmed.rfind("Atoms", 0) == 0){
                insideAtomsSection = true;
                std::getline(input, line);
                continue;
            }
            continue;
        }

        std::istringstream stream(trimmed);
        std::vector<double> fields;
        double field = 0;
        while(stream >> field){
            fields.push_back(field);
        }
        if(fields.size() < 5){
            if(!positions.empty()){
                break;
            }
            continue;
        }

        const int type = static_cast<int>(fields[1]);
        const double x = fields[fields.size() - 3];
        const double y = fields[fields.size() - 2];
        const double z = fields[fields.size() - 1];
        species.push_back(type);
        positions.push_back({x, y, z});
    }

    if(!haveX || !haveY || !haveZ){
        reference.message = "reference file missing box bounds";
        return reference;
    }

    const double lengthA = boxXHigh - boxXLow;
    const double lengthB = boxYHigh - boxYLow;
    const double lengthC = boxZHigh - boxZLow;
    if(lengthA < 1e-6 || lengthB < 1e-6 || lengthC < 1e-6 || positions.empty()){
        reference.message = "degenerate reference cell";
        return reference;
    }

    reference.cellLengthA = lengthA;
    reference.cellLengthB = lengthB;
    reference.cellLengthC = lengthC;

    const double cellLengths[3] = {lengthA, lengthB, lengthC};
    const int siteCount = static_cast<int>(positions.size());
    const int imageRange = static_cast<int>(std::ceil(cutoff / std::min({lengthA, lengthB, lengthC}))) + 1;

    reference.sites.resize(static_cast<std::size_t>(siteCount));
    for(int siteIndex = 0; siteIndex < siteCount; ++siteIndex){
        BasisSite& site = reference.sites[static_cast<std::size_t>(siteIndex)];
        site.species = species[static_cast<std::size_t>(siteIndex)];
        site.fractional = Vector3(
            (positions[siteIndex][0] - boxXLow) / lengthA,
            (positions[siteIndex][1] - boxYLow) / lengthB,
            (positions[siteIndex][2] - boxZLow) / lengthC);

        for(int otherIndex = 0; otherIndex < siteCount; ++otherIndex){
            for(int imageX = -imageRange; imageX <= imageRange; ++imageX){
                for(int imageY = -imageRange; imageY <= imageRange; ++imageY){
                    for(int imageZ = -imageRange; imageZ <= imageRange; ++imageZ){
                        if(siteIndex == otherIndex && imageX == 0 && imageY == 0 && imageZ == 0){
                            continue;
                        }
                        const Vector3 delta(
                            positions[otherIndex][0] + imageX * cellLengths[0] - positions[siteIndex][0],
                            positions[otherIndex][1] + imageY * cellLengths[1] - positions[siteIndex][1],
                            positions[otherIndex][2] + imageZ * cellLengths[2] - positions[siteIndex][2]);
                        const double distance = delta.length();
                        if(distance > 0.1 && distance < cutoff){
                            site.shell.emplace_back(species[static_cast<std::size_t>(otherIndex)], delta);
                        }
                    }
                }
            }
        }
    }

    reference.ok = true;
    return reference;
}

Matrix3 orthonormalFrame(const Vector3& firstVector, const Vector3& secondVector){
    const Vector3 axis1 = firstVector / firstVector.length();
    Vector3 axis2 = secondVector - axis1 * secondVector.dot(axis1);
    axis2 /= axis2.length();
    const Vector3 axis3 = axis1.cross(axis2);
    return Matrix3(axis1, axis2, axis3);
}

std::vector<Vector3> buildIdealNeighborVectors(const CationReference& cationReference){
    constexpr double neighborCutoff = 6.15;
    const double lengthA = cationReference.cellLengthA;
    const double lengthB = cationReference.cellLengthB;
    const double lengthC = cationReference.cellLengthC;

    std::vector<Vector3> replicatedPositions;
    for(const Vector3& fractional : cationReference.idealFractional){
        for(int imageX = -2; imageX <= 2; ++imageX){
            for(int imageY = -2; imageY <= 2; ++imageY){
                for(int imageZ = -2; imageZ <= 2; ++imageZ){
                    replicatedPositions.emplace_back(
                        (fractional.x() + imageX) * lengthA,
                        (fractional.y() + imageY) * lengthB,
                        (fractional.z() + imageZ) * lengthC);
                }
            }
        }
    }

    std::vector<Vector3> idealVectors;
    for(const Vector3& fractional : cationReference.idealFractional){
        const Vector3 center(fractional.x() * lengthA, fractional.y() * lengthB, fractional.z() * lengthC);
        for(const Vector3& replicated : replicatedPositions){
            const Vector3 delta = replicated - center;
            const double distance = delta.length();
            if(distance < 0.1 || distance > neighborCutoff){
                continue;
            }
            bool duplicate = false;
            for(const Vector3& existing : idealVectors){
                if((existing - delta).length() < 0.15){
                    duplicate = true;
                    break;
                }
            }
            if(!duplicate){
                idealVectors.push_back(delta);
            }
        }
    }
    return idealVectors;
}

GrainFrame recoverGrainFrame(
    const Point3* positions,
    const int* types,
    std::size_t atomCount,
    const SimulationCell& cell,
    const CationReference& cationReference)
{
    GrainFrame result;
    constexpr int kNearestNeighbors = 16;
    constexpr double bondCutoff = 6.2;

    std::vector<Point3> cationPositions;
    for(std::size_t atomIndex = 0; atomIndex < atomCount; ++atomIndex){
        if(types[atomIndex] == cationReference.targetSpecies){
            cationPositions.push_back(positions[atomIndex]);
        }
    }
    const int cationCount = static_cast<int>(cationPositions.size());
    result.cationAtomCount = cationCount;
    if(cationCount < 50 || cationReference.idealFractional.empty()){
        result.message = "too few cation atoms or no ideal reference";
        return result;
    }

    const std::vector<Vector3> idealVectors = buildIdealNeighborVectors(cationReference);
    if(idealVectors.size() < 4){
        result.message = "ideal reference produced too few neighbour vectors";
        return result;
    }

    Point3 boundsLow = cationPositions[0];
    Point3 boundsHigh = cationPositions[0];
    for(const Point3& position : cationPositions){
        for(int axis = 0; axis < 3; ++axis){
            boundsLow[axis] = std::min(boundsLow[axis], position[axis]);
            boundsHigh[axis] = std::max(boundsHigh[axis], position[axis]);
        }
    }

    const double gridCellSize = cationReference.physicalBondCutoff;
    const int gridCountX = std::max(1, static_cast<int>((boundsHigh.x() - boundsLow.x()) / gridCellSize) + 1);
    const int gridCountY = std::max(1, static_cast<int>((boundsHigh.y() - boundsLow.y()) / gridCellSize) + 1);
    const int gridCountZ = std::max(1, static_cast<int>((boundsHigh.z() - boundsLow.z()) / gridCellSize) + 1);

    const auto gridCellOf = [&](const Point3& position){
        return std::array<int, 3>{
            std::min(gridCountX - 1, std::max(0, static_cast<int>((position.x() - boundsLow.x()) / gridCellSize))),
            std::min(gridCountY - 1, std::max(0, static_cast<int>((position.y() - boundsLow.y()) / gridCellSize))),
            std::min(gridCountZ - 1, std::max(0, static_cast<int>((position.z() - boundsLow.z()) / gridCellSize)))};
    };
    const auto gridKey = [&](int cellX, int cellY, int cellZ){
        return (static_cast<long long>(cellX) * gridCountY + cellY) * gridCountZ + cellZ;
    };

    std::unordered_map<long long, std::vector<int>> grid;
    for(int cationIndex = 0; cationIndex < cationCount; ++cationIndex){
        const std::array<int, 3> gridCell = gridCellOf(cationPositions[cationIndex]);
        grid[gridKey(gridCell[0], gridCell[1], gridCell[2])].push_back(cationIndex);
    }

    std::vector<std::array<int, kNearestNeighbors>> neighborIndices(static_cast<std::size_t>(cationCount));
    std::vector<std::array<double, kNearestNeighbors>> neighborDistances(static_cast<std::size_t>(cationCount));
    for(int cationIndex = 0; cationIndex < cationCount; ++cationIndex){
        const std::array<int, 3> gridCell = gridCellOf(cationPositions[cationIndex]);
        std::vector<std::pair<double, int>> candidates;
        for(int offsetX = -1; offsetX <= 1; ++offsetX){
            for(int offsetY = -1; offsetY <= 1; ++offsetY){
                for(int offsetZ = -1; offsetZ <= 1; ++offsetZ){
                    const int cellX = gridCell[0] + offsetX;
                    const int cellY = gridCell[1] + offsetY;
                    const int cellZ = gridCell[2] + offsetZ;
                    if(cellX < 0 || cellY < 0 || cellZ < 0 ||
                       cellX >= gridCountX || cellY >= gridCountY || cellZ >= gridCountZ){
                        continue;
                    }
                    const auto bucket = grid.find(gridKey(cellX, cellY, cellZ));
                    if(bucket == grid.end()){
                        continue;
                    }
                    for(int otherIndex : bucket->second){
                        if(otherIndex == cationIndex){
                            continue;
                        }
                        const double distance =
                            cell.wrapVector(cationPositions[otherIndex] - cationPositions[cationIndex]).length();
                        candidates.emplace_back(distance, otherIndex);
                    }
                }
            }
        }
        std::sort(candidates.begin(), candidates.end());
        for(int slot = 0; slot < kNearestNeighbors; ++slot){
            if(slot < static_cast<int>(candidates.size())){
                neighborDistances[static_cast<std::size_t>(cationIndex)][slot] = candidates[static_cast<std::size_t>(slot)].first;
                neighborIndices[static_cast<std::size_t>(cationIndex)][slot] = candidates[static_cast<std::size_t>(slot)].second;
            }else{
                neighborDistances[static_cast<std::size_t>(cationIndex)][slot] = 1e9;
                neighborIndices[static_cast<std::size_t>(cationIndex)][slot] = -1;
            }
        }
    }

    const AffineTransformation& cellMatrix = cell.matrix();
    const double cellSpanX = cellMatrix.column(0).length();
    const double cellSpanY = cellMatrix.column(1).length();
    const double midpointZ = (boundsLow.z() + boundsHigh.z()) * 0.5;

    std::vector<int> bulkIndices;
    for(int cationIndex = 0; cationIndex < cationCount; ++cationIndex){
        const Point3& position = cationPositions[cationIndex];
        if(std::abs(position.z() - midpointZ) < 0.45 * (boundsHigh.z() - boundsLow.z()) &&
           position.x() > boundsLow.x() + 12 && position.x() < boundsLow.x() + cellSpanX - 12 &&
           position.y() > boundsLow.y() + 12 && position.y() < boundsLow.y() + cellSpanY - 12){
            bulkIndices.push_back(cationIndex);
        }
    }
    if(bulkIndices.size() < 20){
        bulkIndices.clear();
        for(int cationIndex = 0; cationIndex < cationCount; ++cationIndex){
            bulkIndices.push_back(cationIndex);
        }
    }

    std::vector<Vector3> bulkBonds;
    for(int cationIndex : bulkIndices){
        for(int slot = 1; slot < kNearestNeighbors; ++slot){
            const int neighbor = neighborIndices[static_cast<std::size_t>(cationIndex)][slot];
            if(neighbor < 0){
                continue;
            }
            const Vector3 bond = cell.wrapVector(cationPositions[neighbor] - cationPositions[cationIndex]);
            if(bond.length() < bondCutoff){
                bulkBonds.push_back(bond);
            }
        }
    }

    const auto refineAndScore = [&](Matrix3 rotation) -> std::pair<Matrix3, double>{
        for(int iteration = 0; iteration < 6; ++iteration){
            double correlation[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
            int matchedCount = 0;
            for(std::size_t bondIndex = 0; bondIndex < bulkBonds.size(); bondIndex += 7){
                const Vector3 measured = bulkBonds[bondIndex];
                const Vector3 inCrystalFrame = rotation.transposed() * measured;
                double bestDistance = 1e9;
                Vector3 bestIdeal = idealVectors[0];
                for(const Vector3& ideal : idealVectors){
                    const double distance = (ideal - inCrystalFrame).length();
                    if(distance < bestDistance){
                        bestDistance = distance;
                        bestIdeal = ideal;
                    }
                }
                if(bestDistance > 1.0){
                    continue;
                }
                correlation[0] += bestIdeal.x() * measured.x();
                correlation[1] += bestIdeal.x() * measured.y();
                correlation[2] += bestIdeal.x() * measured.z();
                correlation[3] += bestIdeal.y() * measured.x();
                correlation[4] += bestIdeal.y() * measured.y();
                correlation[5] += bestIdeal.y() * measured.z();
                correlation[6] += bestIdeal.z() * measured.x();
                correlation[7] += bestIdeal.z() * measured.y();
                correlation[8] += bestIdeal.z() * measured.z();
                ++matchedCount;
            }
            if(matchedCount < 10){
                break;
            }
            double rotationPart[9];
            double stretchPart[9];
            ptm::polar_decomposition_3x3(correlation, true, rotationPart, stretchPart);
            rotation = Matrix3(
                Vector3(rotationPart[0], rotationPart[3], rotationPart[6]),
                Vector3(rotationPart[1], rotationPart[4], rotationPart[7]),
                Vector3(rotationPart[2], rotationPart[5], rotationPart[8]));
        }

        std::vector<double> residuals;
        for(std::size_t bondIndex = 0; bondIndex < bulkBonds.size(); bondIndex += 11){
            const Vector3 inCrystalFrame = rotation.transposed() * bulkBonds[bondIndex];
            double bestDistance = 1e9;
            for(const Vector3& ideal : idealVectors){
                bestDistance = std::min(bestDistance, (ideal - inCrystalFrame).length());
            }
            residuals.push_back(bestDistance);
        }
        return {rotation, medianOf(residuals)};
    };

    double bestResidual = 1e9;
    int seededCount = 0;
    const std::size_t seedStride = std::max<std::size_t>(1, bulkIndices.size() / 6);
    for(std::size_t bulkPosition = 0; bulkPosition < bulkIndices.size() && seededCount < 6; bulkPosition += seedStride){
        const int cationIndex = bulkIndices[bulkPosition];

        std::vector<Vector3> measuredBonds;
        for(int slot = 1; slot < kNearestNeighbors; ++slot){
            const int neighbor = neighborIndices[static_cast<std::size_t>(cationIndex)][slot];
            if(neighbor < 0 || neighborDistances[static_cast<std::size_t>(cationIndex)][slot] > bondCutoff){
                continue;
            }
            measuredBonds.push_back(cell.wrapVector(cationPositions[neighbor] - cationPositions[cationIndex]));
        }
        if(measuredBonds.size() < 6){
            continue;
        }
        ++seededCount;

        const Vector3 measuredFirst = measuredBonds[0];
        Vector3 measuredSecond(0, 0, 0);
        bool haveSecond = false;
        for(std::size_t bondIndex = 1; bondIndex < measuredBonds.size(); ++bondIndex){
            const Vector3 unit = measuredBonds[bondIndex] / measuredBonds[bondIndex].length();
            if(std::abs(unit.dot(measuredFirst / measuredFirst.length())) < 0.85){
                measuredSecond = measuredBonds[bondIndex];
                haveSecond = true;
                break;
            }
        }
        if(!haveSecond){
            continue;
        }

        const Matrix3 measuredFrame = orthonormalFrame(measuredFirst, measuredSecond);
        const double measuredLengthFirst = measuredFirst.length();
        const double measuredLengthSecond = measuredSecond.length();
        const double measuredAngle = std::acos(std::max(-1.0, std::min(1.0,
            measuredFirst.dot(measuredSecond) / (measuredLengthFirst * measuredLengthSecond))));

        for(std::size_t idealFirst = 0; idealFirst < idealVectors.size(); ++idealFirst){
            if(std::abs(idealVectors[idealFirst].length() - measuredLengthFirst) > 0.25){
                continue;
            }
            for(std::size_t idealSecond = 0; idealSecond < idealVectors.size(); ++idealSecond){
                if(idealSecond == idealFirst){
                    continue;
                }
                if(std::abs(idealVectors[idealSecond].length() - measuredLengthSecond) > 0.25){
                    continue;
                }
                const double idealAngle = std::acos(std::max(-1.0, std::min(1.0,
                    idealVectors[idealFirst].dot(idealVectors[idealSecond]) /
                    (idealVectors[idealFirst].length() * idealVectors[idealSecond].length()))));
                if(std::abs(idealAngle - measuredAngle) > 0.12){
                    continue;
                }
                const Matrix3 seedRotation =
                    measuredFrame * orthonormalFrame(idealVectors[idealFirst], idealVectors[idealSecond]).transposed();
                const auto [refinedRotation, residual] = refineAndScore(seedRotation);
                if(residual < bestResidual){
                    bestResidual = residual;
                    result.rotation = refinedRotation;
                }
            }
        }
        if(bestResidual < 0.45){
            break;
        }
    }

    if(bestResidual > 0.9){
        result.message = "grain-basis recovery failed (residual " + std::to_string(bestResidual) + " A)";
        return result;
    }
    result.bulkSnapResidual = bestResidual;
    result.ok = true;
    return result;
}

AllAtomNeighbors buildAllAtomNeighbors(
    const Point3* positions,
    std::size_t atomCount,
    const SimulationCell& cell,
    double cutoff)
{
    AllAtomNeighbors result;
    const int count = static_cast<int>(atomCount);
    result.offsets.assign(static_cast<std::size_t>(count) + 1, 0);

    Point3 boundsLow = positions[0];
    Point3 boundsHigh = positions[0];
    for(int atomIndex = 0; atomIndex < count; ++atomIndex){
        for(int axis = 0; axis < 3; ++axis){
            boundsLow[axis] = std::min(boundsLow[axis], positions[atomIndex][axis]);
            boundsHigh[axis] = std::max(boundsHigh[axis], positions[atomIndex][axis]);
        }
    }

    const double gridCellSize = cutoff;
    const int gridCountX = std::max(1, static_cast<int>((boundsHigh.x() - boundsLow.x()) / gridCellSize) + 1);
    const int gridCountY = std::max(1, static_cast<int>((boundsHigh.y() - boundsLow.y()) / gridCellSize) + 1);
    const int gridCountZ = std::max(1, static_cast<int>((boundsHigh.z() - boundsLow.z()) / gridCellSize) + 1);

    const auto gridCellOf = [&](const Point3& position){
        return std::array<int, 3>{
            std::min(gridCountX - 1, std::max(0, static_cast<int>((position.x() - boundsLow.x()) / gridCellSize))),
            std::min(gridCountY - 1, std::max(0, static_cast<int>((position.y() - boundsLow.y()) / gridCellSize))),
            std::min(gridCountZ - 1, std::max(0, static_cast<int>((position.z() - boundsLow.z()) / gridCellSize)))};
    };
    const auto gridKey = [&](int cellX, int cellY, int cellZ){
        return (static_cast<long long>(cellX) * gridCountY + cellY) * gridCountZ + cellZ;
    };

    std::unordered_map<long long, std::vector<int>> grid;
    for(int atomIndex = 0; atomIndex < count; ++atomIndex){
        const std::array<int, 3> gridCell = gridCellOf(positions[atomIndex]);
        grid[gridKey(gridCell[0], gridCell[1], gridCell[2])].push_back(atomIndex);
    }

    std::vector<std::vector<std::pair<int, Vector3>>> perAtomNeighbors(static_cast<std::size_t>(count));
    const double cutoffSquared = cutoff * cutoff;
    for(int atomIndex = 0; atomIndex < count; ++atomIndex){
        const std::array<int, 3> gridCell = gridCellOf(positions[atomIndex]);
        std::vector<std::pair<double, std::pair<int, Vector3>>> candidates;
        for(int offsetX = -1; offsetX <= 1; ++offsetX){
            for(int offsetY = -1; offsetY <= 1; ++offsetY){
                for(int offsetZ = -1; offsetZ <= 1; ++offsetZ){
                    const int cellX = gridCell[0] + offsetX;
                    const int cellY = gridCell[1] + offsetY;
                    const int cellZ = gridCell[2] + offsetZ;
                    if(cellX < 0 || cellY < 0 || cellZ < 0 ||
                       cellX >= gridCountX || cellY >= gridCountY || cellZ >= gridCountZ){
                        continue;
                    }
                    const auto bucket = grid.find(gridKey(cellX, cellY, cellZ));
                    if(bucket == grid.end()){
                        continue;
                    }
                    for(int otherIndex : bucket->second){
                        if(otherIndex == atomIndex){
                            continue;
                        }
                        const Vector3 delta = cell.wrapVector(positions[otherIndex] - positions[atomIndex]);
                        const double distanceSquared = delta.squaredLength();
                        if(distanceSquared < cutoffSquared && distanceSquared > 0.01){
                            candidates.emplace_back(distanceSquared, std::make_pair(otherIndex, delta));
                        }
                    }
                }
            }
        }
        std::sort(candidates.begin(), candidates.end(), [](const auto& left, const auto& right){
            return left.first < right.first;
        });
        const int keepCount = std::min(static_cast<int>(candidates.size()), kNeighborSlots);
        for(int slot = 0; slot < keepCount; ++slot){
            perAtomNeighbors[static_cast<std::size_t>(atomIndex)].push_back(candidates[static_cast<std::size_t>(slot)].second);
        }
    }

    for(int atomIndex = 0; atomIndex < count; ++atomIndex){
        result.offsets[static_cast<std::size_t>(atomIndex) + 1] =
            result.offsets[static_cast<std::size_t>(atomIndex)] +
            static_cast<int>(perAtomNeighbors[static_cast<std::size_t>(atomIndex)].size());
    }
    result.indices.reserve(static_cast<std::size_t>(result.offsets[static_cast<std::size_t>(count)]));
    result.deltas.reserve(result.indices.capacity());
    for(int atomIndex = 0; atomIndex < count; ++atomIndex){
        for(const auto& [neighbor, delta] : perAtomNeighbors[static_cast<std::size_t>(atomIndex)]){
            result.indices.push_back(neighbor);
            result.deltas.push_back(delta);
        }
    }
    return result;
}

} // namespace

FullCrystalContextResult buildFullCrystalContext(
    const LammpsParser::Frame& frame,
    StructureContext& context,
    StructureAnalysis& analysis,
    const FullCrystalContextParams& params)
{
    FullCrystalContextResult result;
    const std::size_t atomCount = static_cast<std::size_t>(frame.natoms);
    result.atomCount = static_cast<int>(atomCount);

    double referenceCellMin = 0;
    {
        const PerfectReference cellOnly = buildPerfectReference(params.referenceFile, 0.1);
        if(!cellOnly.ok){
            result.message = "perfect reference: " + cellOnly.message;
            return result;
        }
        referenceCellMin = std::min({cellOnly.cellLengthA, cellOnly.cellLengthB, cellOnly.cellLengthC});
    }

    const double referenceReach = std::max(params.bondCutoff, 0.9 * referenceCellMin);
    const PerfectReference reference = buildPerfectReference(params.referenceFile, referenceReach);
    if(!reference.ok){
        result.message = "perfect reference: " + reference.message;
        return result;
    }
    spdlog::info("FullCrystalContext: reference cell=({:.3f},{:.3f},{:.3f}) basis sites={} refReach={:.2f}",
                 reference.cellLengthA, reference.cellLengthB, reference.cellLengthC,
                 reference.sites.size(), referenceReach);

    {
        const double maxCellLength = std::max({reference.cellLengthA, reference.cellLengthB, reference.cellLengthC});
        result.metricRescaleX = reference.cellLengthA / maxCellLength;
        result.metricRescaleY = reference.cellLengthB / maxCellLength;
        result.metricRescaleZ = reference.cellLengthC / maxCellLength;
    }

    CationReference cationReference;
    cationReference.cellLengthA = reference.cellLengthA;
    cationReference.cellLengthB = reference.cellLengthB;
    cationReference.cellLengthC = reference.cellLengthC;
    cationReference.targetSpecies = params.cationSpecies;
    for(const BasisSite& site : reference.sites){
        if(site.species != params.cationSpecies){
            continue;
        }
        Vector3 fractional = site.fractional;
        for(int axis = 0; axis < 3; ++axis){
            fractional[axis] -= std::floor(fractional[axis]);
        }
        bool duplicate = false;
        for(const Vector3& existing : cationReference.idealFractional){
            if((existing - fractional).length() < 0.02){
                duplicate = true;
                break;
            }
        }
        if(!duplicate){
            cationReference.idealFractional.push_back(fractional);
        }
    }

    const GrainFrame grainFrame = recoverGrainFrame(
        frame.positions.data(), frame.types.data(), atomCount, frame.simulationCell, cationReference);
    result.cationAtoms = grainFrame.cationAtomCount;
    if(!grainFrame.ok){
        result.message = "grain frame: " + grainFrame.message;
        return result;
    }
    result.grainSnapResidual = grainFrame.bulkSnapResidual;
    const Matrix3 grainRotation = grainFrame.rotation;
    const Matrix3 grainRotationTransposed = grainRotation.transposed();
    spdlog::info("FullCrystalContext: grain frame recovered (residual {:.3f} A, {} cations)",
                 grainFrame.bulkSnapResidual, grainFrame.cationAtomCount);

    const double rho = ambiguityRho(reference.cellLengthA, reference.cellLengthB, reference.cellLengthC);

    const auto assignForCutoff = [&](
        double cutoff,
        AllAtomNeighbors& neighbors,
        std::vector<Vector3>& overrides,
        FullCrystalContextResult& localResult) -> double
    {
        neighbors = buildAllAtomNeighbors(frame.positions.data(), atomCount, frame.simulationCell, cutoff);

        std::unordered_map<int, std::vector<int>> sitesBySpecies;
        for(int siteIndex = 0; siteIndex < static_cast<int>(reference.sites.size()); ++siteIndex){
            sitesBySpecies[reference.sites[static_cast<std::size_t>(siteIndex)].species].push_back(siteIndex);
        }

        overrides.assign(atomCount * static_cast<std::size_t>(kNeighborSlots), Vector3::Zero());
        localResult.minNeighborCount = std::numeric_limits<int>::max();
        localResult.maxNeighborCount = 0;
        localResult.atomsWithZeroNeighbors = 0;
        localResult.edgesTotal = 0;
        localResult.edgesWithIdealVector = 0;
        localResult.atomsBySpeciesUnassigned = 0;

        const double cellLengths[3] = {reference.cellLengthA, reference.cellLengthB, reference.cellLengthC};

        std::vector<int> basisSiteOfAtom(atomCount, -1);
        for(std::size_t atomIndex = 0; atomIndex < atomCount; ++atomIndex){
            const int start = neighbors.offsets[atomIndex];
            const int end = neighbors.offsets[atomIndex + 1];
            const int neighborCount = end - start;
            localResult.minNeighborCount = std::min(localResult.minNeighborCount, neighborCount);
            localResult.maxNeighborCount = std::max(localResult.maxNeighborCount, neighborCount);
            if(neighborCount == 0){
                ++localResult.atomsWithZeroNeighbors;
                continue;
            }

            const int atomSpecies = frame.types[atomIndex];
            std::vector<Vector3> crystalFrameDeltas(static_cast<std::size_t>(neighborCount));
            for(int slot = 0; slot < neighborCount; ++slot){
                crystalFrameDeltas[static_cast<std::size_t>(slot)] =
                    grainRotationTransposed * neighbors.deltas[static_cast<std::size_t>(start + slot)];
            }

            int bestSite = -1;
            double bestScore = std::numeric_limits<double>::max();
            const auto speciesSites = sitesBySpecies.find(atomSpecies);
            if(speciesSites != sitesBySpecies.end()){
                for(int siteIndex : speciesSites->second){
                    const BasisSite& site = reference.sites[static_cast<std::size_t>(siteIndex)];
                    double score = 0;
                    for(int slot = 0; slot < neighborCount; ++slot){
                        const int neighborSpecies =
                            frame.types[static_cast<std::size_t>(neighbors.indices[static_cast<std::size_t>(start + slot)])];
                        double bestShellDistance = std::numeric_limits<double>::max();
                        for(const auto& [shellSpecies, shellVector] : site.shell){
                            if(shellSpecies != neighborSpecies){
                                continue;
                            }
                            bestShellDistance = std::min(bestShellDistance,
                                (shellVector - crystalFrameDeltas[static_cast<std::size_t>(slot)]).squaredLength());
                        }
                        if(bestShellDistance < std::numeric_limits<double>::max()){
                            score += bestShellDistance;
                        }else{
                            score += cutoff * cutoff;
                        }
                    }
                    if(score < bestScore){
                        bestScore = score;
                        bestSite = siteIndex;
                    }
                }
            }

            localResult.edgesTotal += neighborCount;
            if(bestSite < 0){
                ++localResult.atomsBySpeciesUnassigned;
                continue;
            }
            basisSiteOfAtom[atomIndex] = bestSite;
        }
        if(localResult.minNeighborCount == std::numeric_limits<int>::max()){
            localResult.minNeighborCount = 0;
        }

        std::vector<double> snapResiduals;
        {
            const AffineTransformation& cellMatrix = frame.simulationCell.matrix();
            const double cellSpanX = cellMatrix.column(0).length();
            const double cellSpanY = cellMatrix.column(1).length();
            const double cellSpanZ = cellMatrix.column(2).length();
            Point3 boundsLow = frame.positions[0];
            for(std::size_t atomIndex = 0; atomIndex < atomCount; ++atomIndex){
                for(int axis = 0; axis < 3; ++axis){
                    boundsLow[axis] = std::min(boundsLow[axis], frame.positions[atomIndex][axis]);
                }
            }

            for(std::size_t atomIndex = 0; atomIndex < atomCount; ++atomIndex){
                const int basisSite = basisSiteOfAtom[atomIndex];
                if(basisSite < 0){
                    continue;
                }
                const Vector3& atomFractional = reference.sites[static_cast<std::size_t>(basisSite)].fractional;
                const int start = neighbors.offsets[atomIndex];
                const int end = neighbors.offsets[atomIndex + 1];
                const Point3& atomPosition = frame.positions[atomIndex];
                const bool atomInBulk = !(
                    atomPosition.x() < boundsLow.x() + 0.2 * cellSpanX || atomPosition.x() > boundsLow.x() + 0.8 * cellSpanX ||
                    atomPosition.y() < boundsLow.y() + 0.2 * cellSpanY || atomPosition.y() > boundsLow.y() + 0.8 * cellSpanY ||
                    atomPosition.z() < boundsLow.z() + 0.2 * cellSpanZ || atomPosition.z() > boundsLow.z() + 0.8 * cellSpanZ);

                for(int slot = start; slot < end; ++slot){
                    const int neighbor = neighbors.indices[static_cast<std::size_t>(slot)];
                    const int neighborBasisSite = basisSiteOfAtom[static_cast<std::size_t>(neighbor)];
                    if(neighborBasisSite < 0){
                        continue;
                    }
                    const Vector3& neighborFractional = reference.sites[static_cast<std::size_t>(neighborBasisSite)].fractional;
                    const Vector3 crystalFrameDelta = grainRotationTransposed * neighbors.deltas[static_cast<std::size_t>(slot)];

                    double idealFractionalDelta[3];
                    for(int axis = 0; axis < 3; ++axis){
                        const double measuredFractional = crystalFrameDelta[axis] / cellLengths[axis];
                        const double basisFractional = neighborFractional[axis] - atomFractional[axis];
                        const double imageOffset = std::round(measuredFractional - basisFractional);
                        idealFractionalDelta[axis] = (basisFractional + imageOffset) * cellLengths[axis];
                    }
                    const Vector3 idealLabFrame = grainRotation *
                        Vector3(idealFractionalDelta[0], idealFractionalDelta[1], idealFractionalDelta[2]);

                    overrides[atomIndex * static_cast<std::size_t>(kNeighborSlots) + static_cast<std::size_t>(slot - start)] = idealLabFrame;
                    ++localResult.edgesWithIdealVector;
                    if(atomInBulk && snapResiduals.size() < 20000){
                        snapResiduals.push_back((neighbors.deltas[static_cast<std::size_t>(slot)] - idealLabFrame).length());
                    }
                }
            }
        }

        {
            const auto slotInto = [&](std::size_t fromAtom, int targetAtom) -> int{
                const int start = neighbors.offsets[fromAtom];
                const int end = neighbors.offsets[fromAtom + 1];
                for(int slot = start; slot < end; ++slot){
                    if(neighbors.indices[static_cast<std::size_t>(slot)] == targetAtom){
                        return slot - start;
                    }
                }
                return -1;
            };
            const auto overrideAt = [&](std::size_t atomIndex, int slot) -> Vector3&{
                return overrides[atomIndex * static_cast<std::size_t>(kNeighborSlots) + static_cast<std::size_t>(slot)];
            };

            for(std::size_t atomIndex = 0; atomIndex < atomCount; ++atomIndex){
                const int start = neighbors.offsets[atomIndex];
                const int end = neighbors.offsets[atomIndex + 1];
                for(int slot = start; slot < end; ++slot){
                    const int neighbor = neighbors.indices[static_cast<std::size_t>(slot)];
                    if(neighbor <= static_cast<int>(atomIndex)){
                        continue;
                    }
                    const int slotForward = slot - start;
                    const int slotReverse = slotInto(static_cast<std::size_t>(neighbor), static_cast<int>(atomIndex));
                    if(slotReverse < 0){
                        continue;
                    }
                    Vector3& forwardVector = overrideAt(atomIndex, slotForward);
                    Vector3& reverseVector = overrideAt(static_cast<std::size_t>(neighbor), slotReverse);
                    const bool forwardZero = forwardVector.isZero(1e-9);
                    const bool reverseZero = reverseVector.isZero(1e-9);
                    if(!forwardZero){
                        reverseVector = -forwardVector;
                    }else if(!reverseZero){
                        forwardVector = -reverseVector;
                    }
                }
            }

            localResult.edgesWithIdealVector = 0;
            for(std::size_t atomIndex = 0; atomIndex < atomCount; ++atomIndex){
                const int start = neighbors.offsets[atomIndex];
                const int end = neighbors.offsets[atomIndex + 1];
                for(int slot = start; slot < end; ++slot){
                    if(!overrideAt(atomIndex, slot - start).isZero(1e-9)){
                        ++localResult.edgesWithIdealVector;
                    }
                }
            }
        }

        std::vector<double> loopResiduals;
        {
            const AffineTransformation& cellMatrix = frame.simulationCell.matrix();
            const double cellSpanX = cellMatrix.column(0).length();
            const double cellSpanY = cellMatrix.column(1).length();
            const double cellSpanZ = cellMatrix.column(2).length();
            Point3 boundsLow = frame.positions[0];
            for(std::size_t atomIndex = 0; atomIndex < atomCount; ++atomIndex){
                for(int axis = 0; axis < 3; ++axis){
                    boundsLow[axis] = std::min(boundsLow[axis], frame.positions[atomIndex][axis]);
                }
            }
            const auto overrideOf = [&](std::size_t atomIndex, int slot) -> Vector3{
                return overrides[atomIndex * static_cast<std::size_t>(kNeighborSlots) + static_cast<std::size_t>(slot)];
            };
            const auto findSlot = [&](std::size_t atomIndex, int targetAtom) -> int{
                const int start = neighbors.offsets[atomIndex];
                const int end = neighbors.offsets[atomIndex + 1];
                for(int slot = start; slot < end; ++slot){
                    if(neighbors.indices[static_cast<std::size_t>(slot)] == targetAtom){
                        return slot - start;
                    }
                }
                return -1;
            };

            for(std::size_t atomIndex = 0; atomIndex < atomCount && loopResiduals.size() < 2000; ++atomIndex){
                const Point3& atomPosition = frame.positions[atomIndex];
                if(atomPosition.x() < boundsLow.x() + 0.2 * cellSpanX || atomPosition.x() > boundsLow.x() + 0.8 * cellSpanX){
                    continue;
                }
                if(atomPosition.y() < boundsLow.y() + 0.2 * cellSpanY || atomPosition.y() > boundsLow.y() + 0.8 * cellSpanY){
                    continue;
                }
                if(atomPosition.z() < boundsLow.z() + 0.2 * cellSpanZ || atomPosition.z() > boundsLow.z() + 0.8 * cellSpanZ){
                    continue;
                }

                const int start = neighbors.offsets[atomIndex];
                const int end = neighbors.offsets[atomIndex + 1];
                for(int firstSlot = start; firstSlot < end; ++firstSlot){
                    const int firstNeighbor = neighbors.indices[static_cast<std::size_t>(firstSlot)];
                    for(int secondSlot = firstSlot + 1; secondSlot < end; ++secondSlot){
                        const int secondNeighbor = neighbors.indices[static_cast<std::size_t>(secondSlot)];
                        const int slotFirstToSecond = findSlot(static_cast<std::size_t>(firstNeighbor), secondNeighbor);
                        if(slotFirstToSecond < 0){
                            continue;
                        }
                        const int slotSecondToAtom = findSlot(static_cast<std::size_t>(secondNeighbor), static_cast<int>(atomIndex));
                        if(slotSecondToAtom < 0){
                            continue;
                        }
                        const Vector3 edgeAtomToFirst = overrideOf(atomIndex, firstSlot - start);
                        const Vector3 edgeFirstToSecond = overrideOf(static_cast<std::size_t>(firstNeighbor), slotFirstToSecond);
                        const Vector3 edgeSecondToAtom = overrideOf(static_cast<std::size_t>(secondNeighbor), slotSecondToAtom);
                        if(edgeAtomToFirst.isZero(1e-9) || edgeFirstToSecond.isZero(1e-9) || edgeSecondToAtom.isZero(1e-9)){
                            continue;
                        }
                        loopResiduals.push_back((edgeAtomToFirst + edgeFirstToSecond + edgeSecondToAtom).length());
                        if(loopResiduals.size() >= 2000){
                            break;
                        }
                    }
                    if(loopResiduals.size() >= 2000){
                        break;
                    }
                }
            }
            localResult.bulkLoopsSampled = static_cast<int>(loopResiduals.size());
            localResult.bulkLoopMedianResidual = loopResiduals.empty() ? -1 : medianOf(loopResiduals);
        }

        return snapResiduals.empty() ? std::numeric_limits<double>::max() : percentileOf(snapResiduals, 0.99);
    };

    AllAtomNeighbors neighbors;
    std::vector<Vector3> overrides;
    double chosenCutoff = params.bondCutoff;
    {
        std::vector<double> candidates;
        if(params.bondCutoff > 0.0){
            candidates.push_back(params.bondCutoff);
        }else{
            const double cellMin = std::min({reference.cellLengthA, reference.cellLengthB, reference.cellLengthC});
            const double step = std::max(0.1, cellMin / 24.0);
            const double highestCutoff = 0.85 * cellMin;
            const double lowestCutoff = std::max(step, 2.2);
            for(double cutoff = highestCutoff; cutoff >= lowestCutoff - 1e-9; cutoff -= step){
                candidates.push_back(cutoff);
            }
        }

        bool picked = false;
        FullCrystalContextResult bestLocalResult;
        for(double cutoff : candidates){
            AllAtomNeighbors trialNeighbors;
            std::vector<Vector3> trialOverrides;
            FullCrystalContextResult trialResult;
            const double snapResidualP99 = assignForCutoff(cutoff, trialNeighbors, trialOverrides, trialResult);
            spdlog::info("FullCrystalContext: trial cutoff {:.2f} -> snap p99 {:.3f}A (rho {:.3f}), coord<= {}, idealEdges {}",
                         cutoff, snapResidualP99, rho, trialResult.maxNeighborCount, trialResult.edgesWithIdealVector);
            if(snapResidualP99 < rho || params.bondCutoff > 0.0){
                chosenCutoff = cutoff;
                neighbors = std::move(trialNeighbors);
                overrides = std::move(trialOverrides);
                bestLocalResult = trialResult;
                picked = true;
                break;
            }
        }

        if(!picked){
            const double cutoff = candidates.back();
            FullCrystalContextResult fallbackResult;
            assignForCutoff(cutoff, neighbors, overrides, fallbackResult);
            chosenCutoff = cutoff;
            bestLocalResult = fallbackResult;
            spdlog::warn("FullCrystalContext: no cutoff kept snap p99 < rho {:.3f}A; using tightest {:.2f}A", rho, cutoff);
        }

        result.minNeighborCount = bestLocalResult.minNeighborCount;
        result.maxNeighborCount = bestLocalResult.maxNeighborCount;
        result.atomsWithZeroNeighbors = bestLocalResult.atomsWithZeroNeighbors;
        result.edgesTotal = bestLocalResult.edgesTotal;
        result.edgesWithIdealVector = bestLocalResult.edgesWithIdealVector;
        result.atomsBySpeciesUnassigned = bestLocalResult.atomsBySpeciesUnassigned;
        result.bulkLoopMedianResidual = bestLocalResult.bulkLoopMedianResidual;
        result.bulkLoopsSampled = bestLocalResult.bulkLoopsSampled;
        spdlog::info("FullCrystalContext: SELECTED cutoff {:.2f}A (rho {:.3f}A)", chosenCutoff, rho);
    }
    result.selectedCutoff = chosenCutoff;
    result.ambiguityRho = rho;

    auto offsetsProperty = std::make_shared<ParticleProperty>(atomCount + 1, DataType::Int, 1, 0, true);
    std::copy(neighbors.offsets.begin(), neighbors.offsets.end(), offsetsProperty->dataInt());
    auto indicesProperty = std::make_shared<ParticleProperty>(neighbors.indices.size(), DataType::Int, 1, 0, true);
    if(!neighbors.indices.empty()){
        std::copy(neighbors.indices.begin(), neighbors.indices.end(), indicesProperty->dataInt());
    }
    auto countsProperty = std::make_shared<ParticleProperty>(atomCount, DataType::Int, 1, 0, true);
    for(std::size_t atomIndex = 0; atomIndex < atomCount; ++atomIndex){
        countsProperty->setInt(static_cast<int>(atomIndex), neighbors.offsets[atomIndex + 1] - neighbors.offsets[atomIndex]);
    }

    context.neighborOffsets = offsetsProperty;
    context.neighborIndices = indicesProperty;
    context.neighborCounts = countsProperty;

    result.structureTypesStorage = std::make_shared<ParticleProperty>(atomCount, DataType::Int, 1, 0, true);
    constexpr int kSyntheticStructureType = 1;
    for(std::size_t atomIndex = 0; atomIndex < atomCount; ++atomIndex){
        result.structureTypesStorage->setInt(static_cast<int>(atomIndex), kSyntheticStructureType);
    }
    context.structureTypes = result.structureTypesStorage.get();

    context.atomClusters = std::make_shared<ParticleProperty>(atomCount, DataType::Int, 1, 0, true);
    for(std::size_t atomIndex = 0; atomIndex < atomCount; ++atomIndex){
        context.atomClusters->setInt(static_cast<int>(atomIndex), 1);
    }

    context.maximumNeighborDistance = chosenCutoff;

    analysis.setNeighborLatticeVectorOverrides(std::move(overrides), static_cast<std::size_t>(kNeighborSlots));

    {
        ClusterGraph& clusterGraph = analysis.clusterGraph();
        Cluster* cluster = clusterGraph.createCluster(0, params.topologyName, 1);
        cluster->orientation = grainRotation;
        clusterGraph.createSelfTransition(cluster);
    }

    result.ok = true;
    spdlog::info("FullCrystalContext: atoms={} coord[min={} max={}] zeroNbr={} edges={} idealAssigned={} unassignedSites={} bulkLoopMedian={:.4f}A (n={})",
                 result.atomCount, result.minNeighborCount, result.maxNeighborCount, result.atomsWithZeroNeighbors,
                 result.edgesTotal, result.edgesWithIdealVector, result.atomsBySpeciesUnassigned,
                 result.bulkLoopMedianResidual, result.bulkLoopsSampled);
    return result;
}

}
