#include <channel/io.hpp>

#include <filesystem>
#include <cctype>
#include <fstream>
#include <sstream>
#include <stdexcept>

#include <json-c/json.h>

namespace channel {

namespace fs = std::filesystem;

void ensure_dir(const std::string& path) {
  if (path.empty()) return;
  fs::create_directories(path);
}

static bool parse_two_doubles(const std::string& line, double& a, double& b) {
  std::stringstream ss(line);
  ss >> a >> b;
  return !ss.fail();
}

Profile1D read_profile1d_file(const std::string& path) {
  std::ifstream f(path);
  if (!f) throw std::runtime_error("Cannot open profile file: " + path);

  std::vector<double> z;
  std::vector<double> v;

  std::string line;
  while (std::getline(f, line)) {
    if (line.empty()) continue;
    if (line[0] == '#' || line[0] == ';') continue;
    double a = 0.0, b = 0.0;
    if (!parse_two_doubles(line, a, b)) continue;
    z.push_back(a);
    v.push_back(b);
  }

  if (z.size() < 2) {
    throw std::runtime_error("Profile file has <2 points: " + path);
  }
  return Profile1D(std::move(z), std::move(v));
}

FieldZAlpha read_field_zalpha_file(const std::string& path) {
  std::ifstream f(path);
  if (!f) throw std::runtime_error("Cannot open FieldZAlpha file: " + path);

  std::vector<double> alpha;
  std::vector<double> z;
  std::vector<double> data; // row-major

  std::string line;
  bool have_alpha = false;
  while (std::getline(f, line)) {
    std::string t = line;
    if (t.empty()) continue;

    // alpha header: "# alpha: 0 0.5 1" or "alpha: 0 0.5 1"
    auto trimmed = t;
    // crude trim
    while (!trimmed.empty() && (trimmed.back() == '\r' || trimmed.back() == '\n')) trimmed.pop_back();
    std::string s = trimmed;
    // skip leading spaces
    std::size_t i = 0;
    while (i < s.size() && std::isspace(static_cast<unsigned char>(s[i]))) ++i;
    s = s.substr(i);
    if (s.empty()) continue;

    if (!have_alpha) {
      std::string sh = s;
      if (sh.rfind("#", 0) == 0) sh = sh.substr(1);
      // trim again
      std::size_t j = 0;
      while (j < sh.size() && std::isspace(static_cast<unsigned char>(sh[j]))) ++j;
      sh = sh.substr(j);

      if (sh.rfind("alpha:", 0) == 0) {
        std::stringstream ss(sh.substr(std::string("alpha:").size()));
        double a = 0.0;
        while (ss >> a) alpha.push_back(a);
        if (alpha.size() < 1) throw std::runtime_error("FieldZAlpha: alpha header has no values: " + path);
        have_alpha = true;
        continue;
      }
    }

    if (s[0] == '#' || s[0] == ';') continue;
    // data line: z f(a0) f(a1) ... f(aN-1)
    std::stringstream ss(s);
    double zz = 0.0;
    if (!(ss >> zz)) continue;
    if (!have_alpha) {
      throw std::runtime_error("FieldZAlpha: missing alpha header in file: " + path);
    }
    std::vector<double> row(alpha.size(), 0.0);
    for (std::size_t k = 0; k < alpha.size(); ++k) {
      if (!(ss >> row[k])) {
        throw std::runtime_error("FieldZAlpha: insufficient columns at z=" + std::to_string(zz) + " in " + path);
      }
    }
    z.push_back(zz);
    // append row in row-major order by alpha: we will reorder later
    data.insert(data.end(), row.begin(), row.end());
  }

  if (!have_alpha) throw std::runtime_error("FieldZAlpha: no alpha header in " + path);
  if (z.size() < 2) throw std::runtime_error("FieldZAlpha: <2 z points in " + path);

  // We read data as z-major rows; convert to row-major alpha-major.
  std::size_t Nz = z.size();
  std::size_t Na = alpha.size();
  std::vector<double> data_alpha_major(Nz * Na, 0.0);
  for (std::size_t kz = 0; kz < Nz; ++kz) {
    for (std::size_t ka = 0; ka < Na; ++ka) {
      double v = data[kz * Na + ka];
      data_alpha_major[ka * Nz + kz] = v;
    }
  }

  return FieldZAlpha(std::move(z), std::move(alpha), std::move(data_alpha_major));
}

void write_table(const std::string& path,
                 const std::vector<std::string>& columns,
                 const std::vector<std::vector<double>>& data_columns,
                 const std::string& header_comment) {
  if (columns.size() != data_columns.size()) {
    throw std::runtime_error("write_table: columns and data_columns size mismatch");
  }
  if (columns.empty()) {
    throw std::runtime_error("write_table: empty table");
  }
  std::size_t nrow = data_columns[0].size();
  for (const auto& col : data_columns) {
    if (col.size() != nrow) throw std::runtime_error("write_table: column length mismatch");
  }

  std::ofstream out(path);
  if (!out) throw std::runtime_error("Cannot write table: " + path);

  if (!header_comment.empty()) out << "# " << header_comment << "\n";
  out << "#";
  for (const auto& c : columns) out << " " << c;
  out << "\n";

  for (std::size_t r = 0; r < nrow; ++r) {
    for (std::size_t c = 0; c < columns.size(); ++c) {
      out << data_columns[c][r];
      if (c + 1 < columns.size()) out << " ";
    }
    out << "\n";
  }
}

void write_results_json(const std::string& output_dir, const ResultsIndex& idx) {
  ensure_dir(output_dir);
  std::string path = (fs::path(output_dir) / "results.json").string();

  json_object* root = json_object_new_object();
  json_object_object_add(root, "schema_version", json_object_new_string(idx.schema_version.c_str()));
  json_object_object_add(root, "channel_version", json_object_new_string(idx.channel_version.c_str()));
  json_object_object_add(root, "mode", json_object_new_string(idx.mode.c_str()));
  json_object_object_add(root, "config_used", json_object_new_string(idx.config_used.c_str()));

  json_object* summary = json_object_new_object();
  for (const auto& [k, v] : idx.summary) {
    json_object_object_add(summary, k.c_str(), json_object_new_string(v.c_str()));
  }
  json_object_object_add(root, "summary", summary);

  json_object* datasets = json_object_new_object();
  for (const auto& [name, meta] : idx.datasets) {
    json_object* o = json_object_new_object();
    json_object_object_add(o, "path", json_object_new_string(meta.path.c_str()));
    json_object* cols = json_object_new_array();
    for (const auto& c : meta.columns) json_object_array_add(cols, json_object_new_string(c.c_str()));
    json_object_object_add(o, "columns", cols);
    json_object_object_add(o, "description", json_object_new_string(meta.description.c_str()));
    json_object_object_add(datasets, name.c_str(), o);
  }
  json_object_object_add(root, "datasets", datasets);

  const char* s = json_object_to_json_string_ext(root, JSON_C_TO_STRING_PRETTY);
  std::ofstream out(path);
  if (!out) {
    json_object_put(root);
    throw std::runtime_error("Cannot write results.json: " + path);
  }
  out << s << "\n";
  json_object_put(root);
}

void copy_file(const std::string& src, const std::string& dst) {
  std::ifstream in(src, std::ios::binary);
  if (!in) throw std::runtime_error("copy_file: cannot open src " + src);
  std::ofstream out(dst, std::ios::binary);
  if (!out) throw std::runtime_error("copy_file: cannot open dst " + dst);
  out << in.rdbuf();
}

} // namespace channel
