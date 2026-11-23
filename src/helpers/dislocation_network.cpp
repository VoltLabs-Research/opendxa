#include <volt/core/volt.h>
#include <volt/core/progress.h>
#include <volt/helpers/dislocation_network.h>
#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
#include <numeric>

namespace Volt{

DislocationSegment* DislocationNetwork::createSegment(const ClusterVector& burgersVector){
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

namespace{

// Replaces runs of consecutive points with their mean. A run grows while its point count squared stays below
// linePointInterval times its summed core size, so lines with wide cores keep fewer points. An open line keeps its ends.
std::vector<Point3> coarsenDislocationLine(
	double linePointInterval,
	const std::vector<Point3>& input,
	const std::vector<int>& coreSize,
	bool isClosedLoop,
	bool isInfiniteLine
){
	if(input.size() < 2){
		return {};
	}

	if(linePointInterval <= 0 || input.size() < 4){
		return input;
	}

	// An infinite line that is short against its core collapses to a straight line through its mean point.
	if(isInfiniteLine){
		const int coreSizeSum = std::accumulate(coreSize.cbegin(), coreSize.cend() - 1, 0);
		const int count = static_cast<int>(input.size()) - 1;
		if(coreSizeSum * linePointInterval > count * count){
			Vector3 com = Vector3::Zero();
			for(auto p = input.cbegin(); p != input.cend() - 1; ++p){
				com += (*p - input.front());
			}
			return {input.front() + com / count, input.back() + com / count};
		}
	}

	std::vector<Point3> output;
	if(!isClosedLoop){
		output.push_back(input.front());
	}

	const size_t minNumPoints = isClosedLoop ? 4 : 2;

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

		output.push_back(Point3::Origin() + com / count);
	}

	if(!isClosedLoop){
		output.push_back(input.back());
	}else{
		output.push_back(output.front());
	}

	if(output.size() < minNumPoints){
		return input;
	}
	return output;
}

void smoothDislocationLine(double smoothingLevel, std::vector<Point3>& line, bool isLoop){
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

}

void DislocationNetwork::smoothDislocationLines(double lineSmoothingLevel, double linePointInterval){
	const auto& segmentList = segments();
	ProgressReporter progress("Smoothing dislocation lines", segmentList.size());
	tbb::parallel_for(tbb::blocked_range<size_t>(0, segmentList.size(), 128),
		[&](const tbb::blocked_range<size_t>& r){
		for(size_t i = r.begin(); i < r.end(); ++i){
			DislocationSegment* segment = segmentList[i];
			if(!segment || segment->coreSize.empty()) continue;

			std::vector<Point3> line = coarsenDislocationLine(
				linePointInterval,
				segment->line,
				segment->coreSize,
				segment->isClosedLoop(),
				segment->isInfiniteLine()
			);
			smoothDislocationLine(lineSmoothingLevel, line, segment->isClosedLoop());
			segment->line = std::move(line);

			segment->coreSize.clear();
		}
		progress.advance(r.size());
	});
}

}
