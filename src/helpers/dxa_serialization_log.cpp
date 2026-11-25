#include <volt/helpers/dxa_serialization_log.h>

#include <volt/helpers/dxa_serialization_detail.h>

#include <spdlog/spdlog.h>

#include <string>

namespace Volt{

void DxaSerializationLog::dump(const DxaSerialization::Detail::LineSummary& summary){
    if(summary.skippedSegments > 0){
        spdlog::info("  Skipped {} degenerate segments (fewer than two points)", summary.skippedSegments);
    }
    spdlog::info(
        "  {} dislocation lines{} ({} points), length min {:.2f} / mean {:.2f} / max {:.2f}",
        summary.lineCount, summary.clipped ? " after PBC clipping" : "", summary.totalPoints,
        summary.minLength, summary.meanLength, summary.maxLength
    );
    if(summary.byFamily.empty()) return;

    std::string families;
    for(const auto& [_, tally] : summary.byFamily){
        if(!families.empty()) families += ", ";
        families += spdlog::fmt_lib::format("{} {} ({:.1f}%)", tally.label, tally.segmentCount,
            summary.totalLength > 0 ? 100.0 * tally.totalLength / summary.totalLength : 0.0);
    }
    spdlog::info("  Burgers vector families (lines, % of length): {}", families);
}

}
