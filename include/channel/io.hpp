#pragma once

#include <channel/profile.hpp>

#include <map>
#include <string>
#include <vector>

namespace channel {

struct DatasetMeta {
  std::string path;                  // relative to output_dir
  std::vector<std::string> columns;  // e.g. ["z","psi"]
  std::string description;
};

struct ResultsIndex {
  std::string schema_version = "channel.results.v1";
  std::string channel_version = "0.1.0";
  std::string mode;          // stationary | dynamic | both
  std::string config_used;   // usually "config_used.ini"

  // High-level run summary (stringified numbers to keep it simple).
  std::map<std::string, std::string> summary;

  // Named datasets produced by this run.
  std::map<std::string, DatasetMeta> datasets;
};

void ensure_dir(const std::string& path);

// Read 2-column profile: z value
Profile1D read_profile1d_file(const std::string& path);

// Read FieldZAlpha file format:
//   # alpha: 0 0.5 1
//   z  f(a0) f(a1) f(a2)
//   ...
FieldZAlpha read_field_zalpha_file(const std::string& path);

// Write a whitespace table with optional header lines beginning with '#'.
void write_table(const std::string& path,
                 const std::vector<std::string>& columns,
                 const std::vector<std::vector<double>>& data_columns,
                 const std::string& header_comment = "");

// Write results.json in output_dir.
void write_results_json(const std::string& output_dir, const ResultsIndex& idx);

// Copy a file (best-effort; overwrites).
void copy_file(const std::string& src, const std::string& dst);

} // namespace channel
