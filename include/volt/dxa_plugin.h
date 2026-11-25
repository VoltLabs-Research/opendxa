#pragma once

#include <volt/plugin/plugin_main.h>

#include <string>
#include <vector>

namespace Volt::DxaPlugin {

Plugin::PluginDescriptor buildDescriptor();

// Offers the topology names of the YAMLs in --lattice_dir as the values of --reference_topology.
void refineTopologyDomain(std::vector<Plugin::CliOption>& options, const Plugin::OptsMap& opts);

// The lattice YAML's own name for a known topology. A name no YAML knows comes from the upstream contract and is used
// as given.
std::string canonicalTopologyName(const std::string& name);

}
