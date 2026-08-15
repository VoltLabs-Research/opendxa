#include <volt/core/volt.h>
#include <volt/helpers/dislocation_network.h>
#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>

namespace Volt{

DislocationNetwork::DislocationNetwork(const DislocationNetwork& other): _clusterGraph(other._clusterGraph){
	_segments.reserve(other._segments.size());

	for(const auto *oldSegment : other.segments()){
		assert(oldSegment->replacedWith == nullptr);
		assert(oldSegment->id == static_cast<int>(_segments.size()));

		auto* newSegment = createSegment(oldSegment->burgersVector);
		newSegment->line = oldSegment->line;
		newSegment->coreSize = oldSegment->coreSize;

		assert(newSegment->id == oldSegment->id);
	}

	for(size_t segmentIndex = 0; segmentIndex < other.segments().size(); ++segmentIndex){
		const auto *oldSegment = other.segments()[segmentIndex];
		auto *newSegment = _segments[segmentIndex];

		for(int nodeIndex = 0; nodeIndex < 2; ++nodeIndex){
			const auto *oldNode = oldSegment->nodes[nodeIndex];
			if(oldNode->isDangling()) continue;

			auto *oldSecondNode = oldNode->junctionRing;
			auto *newNode = newSegment->nodes[nodeIndex];
			auto *newSecondNode= _segments[oldSecondNode->segment->id]->nodes[oldSecondNode->isForwardNode() ? 0 : 1];
			newNode->junctionRing = newSecondNode;
		}
	}
}

DislocationSegment* DislocationNetwork::createSegment(const ClusterVector& burgersVector){
	tbb::spin_mutex::scoped_lock lock(_segmentsMutex);
	DislocationNode *forwardNode = _nodePool.construct();
	DislocationNode *backwardNode = _nodePool.construct();

	DislocationSegment *segment = _segmentPool.construct(burgersVector, forwardNode, backwardNode);
	segment->id = _segments.size();
	_segments.push_back(segment);

	return segment;
}

void DislocationNetwork::discardSegment(DislocationSegment* segment){
	assert(segment != nullptr);
	const auto it = std::ranges::find(_segments, segment);
	assert(it != _segments.end());
	_segments.erase(it);
}

void DislocationNetwork::smoothDislocationLines(double lineSmoothingLevel, double linePointInterval){
	const auto& segmentList = segments();
	tbb::parallel_for(tbb::blocked_range<size_t>(0, segmentList.size(), 128),
		[&](const tbb::blocked_range<size_t>& r){
		for(size_t i = r.begin(); i < r.end(); ++i){
			DislocationSegment* segment = segmentList[i];
			if(!segment || segment->coreSize.empty()) continue;

			std::vector<Point3> line;
			std::vector<int> coreSize;
			coarsenDislocationLine(
				linePointInterval,
				segment->line,
				segment->coreSize,
				line,
				coreSize,
				segment->isClosedLoop(),
				segment->isInfiniteLine()
			);
			smoothDislocationLine(lineSmoothingLevel, line, segment->isClosedLoop());
			segment->line = std::move(line);

			segment->coreSize.clear();
		}
	});
}

void DislocationNetwork::coarsenDislocationLine(
    double linePointInterval,
    const std::vector<Point3>& input,
    const std::vector<int>& coreSize,
    std::vector<Point3>& output,
    std::vector<int>& outputCoreSize,
    bool isClosedLoop,
    bool isInfiniteLine
){
	if(input.size() < 2){
		return;
	}

    if(linePointInterval <= 0 || input.size() < 4){
        output = input;
        outputCoreSize = coreSize;
        return;
    }

    if(isInfiniteLine && input.size() >= 3){
        int coreSizeSum = std::accumulate(coreSize.cbegin(), coreSize.cend() - 1, 0);
        int count = input.empty() ? 0 : static_cast<int>(input.size()) - 1;
        if(count > 0 && coreSizeSum * linePointInterval > count * count){
		Vector3 com = Vector3::Zero();
            for(auto p = input.cbegin(); p != input.cend() - 1; ++p){
                com += (*p - input.front());
            }
            if(count != 0){
                output.push_back(input.front() + com / count);
                outputCoreSize.push_back(coreSizeSum / count);
                output.push_back(input.back() + com / count);
                outputCoreSize.push_back(coreSizeSum / count);
            }
            return;
        }
    }

    output.clear();
    outputCoreSize.clear();

    if(!isClosedLoop){
        output.push_back(input.front());
        outputCoreSize.push_back(coreSize.front());
    }

    size_t minNumPoints = isClosedLoop ? 4 : 2;

    auto point_it = isClosedLoop ? input.cbegin() : std::next(input.cbegin());
    auto core_it = isClosedLoop ? coreSize.cbegin() : std::next(coreSize.cbegin());
    auto end_it = isClosedLoop ? input.cend() : std::prev(input.cend());

    while(point_it != end_it){
		Vector3 com = Vector3::Zero();
        int sum = 0;
        int count = 0;

        do{
		com += (*point_it - Point3::Origin());
            sum += *core_it;
            count++;
            ++point_it;
            ++core_it;
        }while(
            point_it != end_it &&
            (static_cast<double>(count) * count) < (linePointInterval * sum)
        );

        if(count > 0){
		output.push_back(Point3::Origin() + com / count);
            outputCoreSize.push_back(sum / count);
        }
    }

    if(!isClosedLoop){
        output.push_back(input.back());
        outputCoreSize.push_back(coreSize.back());
    }else if(!output.empty()){
        output.push_back(output.front());
        outputCoreSize.push_back(outputCoreSize.front());
    }

    if(output.size() < minNumPoints){
        output = input;
        outputCoreSize = coreSize;
    }
}
void DislocationNetwork::smoothDislocationLine(double smoothingLevel, std::vector<Point3>& line, bool isLoop){
	if(smoothingLevel <= 0 || line.size() <= 2){
		return;
	}

	if(isLoop && line.size() <= 4){
		return;
	}

    const double lambda = 0.5;
    const double mu = -0.53; 
    const double prefactors[2] = { lambda, mu };
    
    std::vector<Vector3> laplacians(line.size());

	for(int iteration = 0; static_cast<double>(iteration) < smoothingLevel; ++iteration){
        for(int pass = 0; pass < 2; ++pass){
			
			for(size_t i = 0; i < line.size(); ++i){
                if(!isLoop && (i == 0 || i == line.size() - 1)){
					laplacians[i].setZero();
					continue;
				}

				const Point3& p_prev = (i == 0) ? line[line.size() - 2] : line[i - 1];
				const Point3& p_curr = line[i];
				const Point3& p_next = (i == line.size() - 1) ? line[1] : line[i + 1];

                laplacians[i] = 0.5 * ((p_prev - p_curr) + (p_next - p_curr));
			}

			for(size_t i = 0; i < line.size(); ++i){
				line[i] += prefactors[pass] * laplacians[i];
			}
		}
	}
}
[[nodiscard]] Point3 DislocationSegment::getPointOnLine(double t) const{
	if(line.empty() || isDegenerate()){
				return Point3::Origin();
	}

	t *= calculateLength();
	double sum = 0;

	for(auto i1 = line.begin(); std::next(i1) != line.end(); ++i1){
		const auto i2 = std::next(i1);
		const auto delta = *i2 - *i1;
		const auto len = delta.length();

		if(len != 0 && sum + len >= t){
			return *i1 + ((t - sum) / len) * delta;
		}

		sum += len;
	}

	return line.back();
}

}
