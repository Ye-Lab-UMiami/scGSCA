// ============================================================================
// DataFrame.h
//
// Lightweight columnar data frame for the MetaGSCA pipeline.
// Each column is stored as a variant: either vector<double> (numeric) or
// vector<string> (character / factor).  Row-subsetting, column-subsetting,
// cbind, and rbind are supported.
// ============================================================================
#ifndef METAGSCA_DATAFRAME_H
#define METAGSCA_DATAFRAME_H

#include <string>
#include <vector>
#include <unordered_map>
#include <variant>
#include <stdexcept>
#include <algorithm>
#include <cmath>
#include <cassert>

namespace gsca {

// A column is either numeric or string-valued.
using NumCol = std::vector<double>;
using StrCol = std::vector<std::string>;
using Column = std::variant<NumCol, StrCol>;

// -------------------------------------------------------------------------
// DataFrame
// -------------------------------------------------------------------------
class DataFrame {
public:
    DataFrame() = default;

    // --- dimensions -------------------------------------------------------
    size_t nrow() const { return nrow_; }
    size_t ncol() const { return col_names_.size(); }

    // --- column access by name --------------------------------------------
    bool has_col(const std::string& name) const;
    const Column& col(const std::string& name) const;
    Column&       col(const std::string& name);

    const NumCol& num_col(const std::string& name) const;
    NumCol&       num_col(const std::string& name);

    const StrCol& str_col(const std::string& name) const;
    StrCol&       str_col(const std::string& name);

    bool is_numeric(const std::string& name) const;

    // --- column access by index -------------------------------------------
    const std::string& col_name(size_t idx) const { return col_names_.at(idx); }
    const Column& col_by_idx(size_t idx) const { return columns_.at(idx); }

    // --- column names -----------------------------------------------------
    const std::vector<std::string>& col_names() const { return col_names_; }

    // --- row names (optional, like R rownames) ----------------------------
    std::vector<std::string> row_names;

    // --- add columns ------------------------------------------------------
    void add_num_col(const std::string& name, NumCol data);
    void add_str_col(const std::string& name, StrCol data);
    void add_col(const std::string& name, Column data);

    // --- remove column ----------------------------------------------------
    void remove_col(const std::string& name);
    void rename_col(const std::string& old_name, const std::string& new_name);

    // --- subset rows (returns new DataFrame) ------------------------------
    DataFrame row_subset(const std::vector<size_t>& indices) const;
    DataFrame row_subset(const std::vector<bool>& mask) const;

    // --- subset columns (returns new DataFrame) ---------------------------
    DataFrame col_subset(const std::vector<std::string>& names) const;

    // --- rbind (append rows from another DataFrame) -----------------------
    static DataFrame rbind(const DataFrame& a, const DataFrame& b);

    // --- cbind (prepend a numeric column; used for cid) -------------------
    static DataFrame cbind_num(const std::string& name, const NumCol& data,
                               const DataFrame& df);

    // --- round numeric columns in place -----------------------------------
    void round_cols(size_t from_col_idx);

    // --- replace NA (NaN) with 0 in numeric columns -----------------------
    void replace_nan(size_t from_col_idx, double replacement = 0.0);

    // --- gene column names from gene.start to end -------------------------
    std::vector<std::string> gene_names(size_t gene_start_idx) const;

private:
    size_t nrow_ = 0;
    std::vector<std::string> col_names_;
    std::vector<Column> columns_;
    std::unordered_map<std::string, size_t> name_to_idx_;

    void check_length(size_t len) const;
    size_t idx_of(const std::string& name) const;
};

// -------------------------------------------------------------------------
// Utility: read a numeric column as Eigen::VectorXd (declared in .cpp)
// -------------------------------------------------------------------------

} // namespace gsca

#endif // METAGSCA_DATAFRAME_H
