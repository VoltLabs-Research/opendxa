#pragma once

#include <volt/core/run_log.h>

namespace Volt{

namespace DxaSerialization::Detail{
struct LineSummary;
}

// The log of the dislocation export: what the lines written for a network add up to. Scripts parse these lines, so
// their wording stays as it is.
class DxaSerializationLog : public RunLog{
public:
    using RunLog::dump;

    void dump(const DxaSerialization::Detail::LineSummary& summary);
};

}
