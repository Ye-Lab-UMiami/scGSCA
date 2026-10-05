// ============================================================================
// tsv_io.cpp — Read / write tab-separated files.
// ============================================================================
#include "tsv_io.h"
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <iostream>
#include <cmath>
#include <algorithm>

namespace gsca {

// ---------------------------------------------------------------------------
// Helper: strip leading/trailing whitespace and quotes from a string.
// Also removes carriage return (\r) which may appear in Windows line endings.
// ---------------------------------------------------------------------------
static std::string strip(const std::string& s) {
    std::string out = s;
    // Remove \r
    out.erase(std::remove(out.begin(), out.end(), '\r'), out.end());
    // Trim whitespace
    size_t start = out.find_first_not_of(" \t\n");
    if (start == std::string::npos) return "";
    size_t end = out.find_last_not_of(" \t\n");
    out = out.substr(start, end - start + 1);
    // Remove surrounding quotes
    if (out.size() >= 2 &&
        ((out.front() == '"' && out.back() == '"') ||
         (out.front() == '\'' && out.back() == '\''))) {
        out = out.substr(1, out.size() - 2);
    }
    return out;
}

// ---------------------------------------------------------------------------
// Helper: remove UTF-8 BOM if present at the start of a string.
// ---------------------------------------------------------------------------
static void strip_bom(std::string& s) {
    if (s.size() >= 3 &&
        static_cast<unsigned char>(s[0]) == 0xEF &&
        static_cast<unsigned char>(s[1]) == 0xBB &&
        static_cast<unsigned char>(s[2]) == 0xBF) {
        s = s.substr(3);
    }
}

// ---------------------------------------------------------------------------
// Helper: try to parse every element of a string vector as double.
// Returns true if all elements are valid doubles.
// ---------------------------------------------------------------------------
static bool try_parse_numeric(const std::vector<std::string>& strs,
                              std::vector<double>& out) {
    out.resize(strs.size());
    for (size_t i = 0; i < strs.size(); ++i) {
        if (strs[i].empty() || strs[i] == "NA" || strs[i] == "NaN") {
            out[i] = std::nan("");
            continue;
        }
        try {
            size_t pos = 0;
            out[i] = std::stod(strs[i], &pos);
            if (pos != strs[i].size()) return false;
        } catch (...) {
            return false;
        }
    }
    return true;
}

// ---------------------------------------------------------------------------
// read_tsv
// ---------------------------------------------------------------------------
DataFrame read_tsv(const std::string& path) {
    std::ifstream ifs(path);
    if (!ifs.is_open())
        throw std::runtime_error("read_tsv: cannot open file: " + path);

    // Read header line
    std::string header_line;
    if (!std::getline(ifs, header_line))
        throw std::runtime_error("read_tsv: empty file: " + path);

    // Remove BOM if present
    strip_bom(header_line);

    // Remove trailing \r (Windows line endings)
    header_line.erase(std::remove(header_line.begin(), header_line.end(), '\r'),
                      header_line.end());

    // Parse column names
    std::vector<std::string> col_names;
    {
        std::istringstream hs(header_line);
        std::string tok;
        while (std::getline(hs, tok, '\t'))
            col_names.push_back(strip(tok));
    }

    // Debug: print detected column names
    std::cout << "Detected columns (" << col_names.size() << "): ";
    for (size_t i = 0; i < std::min(col_names.size(), size_t(6)); ++i) {
        if (i > 0) std::cout << ", ";
        std::cout << "'" << col_names[i] << "'";
    }
    if (col_names.size() > 6) std::cout << " ...";
    std::cout << std::endl;

    size_t n_cols = col_names.size();
    // Temporary storage: all columns as strings
    std::vector<std::vector<std::string>> str_cols(n_cols);

    std::string line;
    while (std::getline(ifs, line)) {
        // Remove \r
        line.erase(std::remove(line.begin(), line.end(), '\r'), line.end());
        if (line.empty()) continue;

        std::istringstream ls(line);
        std::string tok;
        size_t c = 0;
        while (std::getline(ls, tok, '\t') && c < n_cols) {
            str_cols[c].push_back(strip(tok));
            ++c;
        }
        // Pad missing columns with empty strings
        while (c < n_cols) {
            str_cols[c].push_back("");
            ++c;
        }
    }

    // Build DataFrame, auto-detecting numeric columns
    DataFrame df;
    for (size_t c = 0; c < n_cols; ++c) {
        std::vector<double> nums;
        if (try_parse_numeric(str_cols[c], nums)) {
            df.add_num_col(col_names[c], std::move(nums));
        } else {
            df.add_str_col(col_names[c], std::move(str_cols[c]));
        }
    }

    // Generate default row names
    df.row_names.resize(df.nrow());
    for (size_t i = 0; i < df.nrow(); ++i)
        df.row_names[i] = std::to_string(i + 1);

    return df;
}

// ---------------------------------------------------------------------------
// write_tsv
// ---------------------------------------------------------------------------
void write_tsv(const DataFrame& df, const std::string& path) {
    std::ofstream ofs(path);
    if (!ofs.is_open())
        throw std::runtime_error("write_tsv: cannot open file: " + path);

    // Header
    for (size_t c = 0; c < df.ncol(); ++c) {
        if (c > 0) ofs << '\t';
        ofs << df.col_name(c);
    }
    ofs << '\n';

    // Data
    for (size_t r = 0; r < df.nrow(); ++r) {
        for (size_t c = 0; c < df.ncol(); ++c) {
            if (c > 0) ofs << '\t';
            std::visit([&](auto& vec) {
                using T = std::decay_t<decltype(vec)>;
                if constexpr (std::is_same_v<T, NumCol>) {
                    double v = vec[r];
                    if (std::isnan(v))
                        ofs << "NA";
                    else
                        ofs << v;
                } else {
                    ofs << vec[r];
                }
            }, df.col_by_idx(c));
        }
        ofs << '\n';
    }
}

} // namespace gsca
