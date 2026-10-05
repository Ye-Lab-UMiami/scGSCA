// ============================================================================
// DataFrame.cpp — implementation of the lightweight columnar data frame.
// ============================================================================
#include "DataFrame.h"
#include <sstream>

namespace gsca {

// --- helpers ---------------------------------------------------------------
void DataFrame::check_length(size_t len) const {
    if (nrow_ > 0 && len != nrow_) {
        throw std::runtime_error(
            "DataFrame: column length (" + std::to_string(len) +
            ") != existing nrow (" + std::to_string(nrow_) + ")");
    }
}

size_t DataFrame::idx_of(const std::string& name) const {
    auto it = name_to_idx_.find(name);
    if (it == name_to_idx_.end())
        throw std::runtime_error("DataFrame: column '" + name + "' not found.");
    return it->second;
}

// --- has_col ---------------------------------------------------------------
bool DataFrame::has_col(const std::string& name) const {
    return name_to_idx_.count(name) > 0;
}

// --- column access ---------------------------------------------------------
const Column& DataFrame::col(const std::string& name) const {
    return columns_.at(idx_of(name));
}
Column& DataFrame::col(const std::string& name) {
    return columns_.at(idx_of(name));
}

const NumCol& DataFrame::num_col(const std::string& name) const {
    return std::get<NumCol>(col(name));
}
NumCol& DataFrame::num_col(const std::string& name) {
    return std::get<NumCol>(col(name));
}

const StrCol& DataFrame::str_col(const std::string& name) const {
    return std::get<StrCol>(col(name));
}
StrCol& DataFrame::str_col(const std::string& name) {
    return std::get<StrCol>(col(name));
}

bool DataFrame::is_numeric(const std::string& name) const {
    return std::holds_alternative<NumCol>(col(name));
}

// --- add columns -----------------------------------------------------------
void DataFrame::add_num_col(const std::string& name, NumCol data) {
    check_length(data.size());
    if (nrow_ == 0) nrow_ = data.size();
    name_to_idx_[name] = columns_.size();
    col_names_.push_back(name);
    columns_.emplace_back(std::move(data));
}

void DataFrame::add_str_col(const std::string& name, StrCol data) {
    check_length(data.size());
    if (nrow_ == 0) nrow_ = data.size();
    name_to_idx_[name] = columns_.size();
    col_names_.push_back(name);
    columns_.emplace_back(std::move(data));
}

void DataFrame::add_col(const std::string& name, Column data) {
    size_t len = std::visit([](auto& v) { return v.size(); }, data);
    check_length(len);
    if (nrow_ == 0) nrow_ = len;
    name_to_idx_[name] = columns_.size();
    col_names_.push_back(name);
    columns_.emplace_back(std::move(data));
}

// --- remove column ---------------------------------------------------------
void DataFrame::remove_col(const std::string& name) {
    size_t idx = idx_of(name);
    columns_.erase(columns_.begin() + static_cast<long>(idx));
    col_names_.erase(col_names_.begin() + static_cast<long>(idx));
    // rebuild index
    name_to_idx_.clear();
    for (size_t i = 0; i < col_names_.size(); ++i)
        name_to_idx_[col_names_[i]] = i;
}

void DataFrame::rename_col(const std::string& old_name, const std::string& new_name) {
    if (old_name == new_name) return;
    size_t idx = idx_of(old_name);
    if (has_col(new_name))
        throw std::runtime_error("DataFrame: column '" + new_name + "' already exists.");
    col_names_[idx] = new_name;
    name_to_idx_.erase(old_name);
    name_to_idx_[new_name] = idx;
}

// --- row subset by indices -------------------------------------------------
DataFrame DataFrame::row_subset(const std::vector<size_t>& indices) const {
    DataFrame out;
    for (size_t c = 0; c < ncol(); ++c) {
        std::visit([&](auto& vec) {
            using T = std::decay_t<decltype(vec)>;
            T sub;
            sub.reserve(indices.size());
            for (auto i : indices) sub.push_back(vec.at(i));
            out.add_col(col_names_[c], Column(std::move(sub)));
        }, columns_[c]);
    }
    // subset row_names if present
    if (!row_names.empty()) {
        out.row_names.reserve(indices.size());
        for (auto i : indices) out.row_names.push_back(row_names.at(i));
    }
    return out;
}

// --- row subset by boolean mask --------------------------------------------
DataFrame DataFrame::row_subset(const std::vector<bool>& mask) const {
    std::vector<size_t> idx;
    for (size_t i = 0; i < mask.size(); ++i)
        if (mask[i]) idx.push_back(i);
    return row_subset(idx);
}

// --- column subset ---------------------------------------------------------
DataFrame DataFrame::col_subset(const std::vector<std::string>& names) const {
    DataFrame out;
    for (auto& n : names) {
        out.add_col(n, col(n));
    }
    out.row_names = row_names;
    return out;
}

// --- rbind -----------------------------------------------------------------
DataFrame DataFrame::rbind(const DataFrame& a, const DataFrame& b) {
    if (a.ncol() != b.ncol())
        throw std::runtime_error("rbind: ncol mismatch.");
    DataFrame out;
    for (size_t c = 0; c < a.ncol(); ++c) {
        const auto& name = a.col_name(c);
        std::visit([&](auto& va) {
            using T = std::decay_t<decltype(va)>;
            auto& vb = std::get<T>(b.col_by_idx(c));
            T combined;
            combined.reserve(va.size() + vb.size());
            combined.insert(combined.end(), va.begin(), va.end());
            combined.insert(combined.end(), vb.begin(), vb.end());
            out.add_col(name, Column(std::move(combined)));
        }, a.col_by_idx(c));
    }
    // merge row names
    out.row_names.reserve(a.nrow() + b.nrow());
    if (!a.row_names.empty()) {
        out.row_names.insert(out.row_names.end(),
                             a.row_names.begin(), a.row_names.end());
    } else {
        for (size_t i = 0; i < a.nrow(); ++i)
            out.row_names.push_back(std::to_string(i));
    }
    if (!b.row_names.empty()) {
        out.row_names.insert(out.row_names.end(),
                             b.row_names.begin(), b.row_names.end());
    } else {
        for (size_t i = 0; i < b.nrow(); ++i)
            out.row_names.push_back(std::to_string(a.nrow() + i));
    }
    return out;
}

// --- cbind_num (prepend a numeric column) ----------------------------------
DataFrame DataFrame::cbind_num(const std::string& name, const NumCol& data,
                               const DataFrame& df) {
    DataFrame out;
    out.add_num_col(name, data);
    for (size_t c = 0; c < df.ncol(); ++c)
        out.add_col(df.col_name(c), df.col_by_idx(c));
    out.row_names = df.row_names;
    return out;
}

// --- round numeric columns -------------------------------------------------
void DataFrame::round_cols(size_t from_col_idx) {
    for (size_t c = from_col_idx; c < ncol(); ++c) {
        if (std::holds_alternative<NumCol>(columns_[c])) {
            auto& v = std::get<NumCol>(columns_[c]);
            for (auto& x : v) x = std::round(x);
        }
    }
}

// --- replace NaN -----------------------------------------------------------
void DataFrame::replace_nan(size_t from_col_idx, double replacement) {
    for (size_t c = from_col_idx; c < ncol(); ++c) {
        if (std::holds_alternative<NumCol>(columns_[c])) {
            auto& v = std::get<NumCol>(columns_[c]);
            for (auto& x : v)
                if (std::isnan(x)) x = replacement;
        }
    }
}

// --- gene_names ------------------------------------------------------------
std::vector<std::string> DataFrame::gene_names(size_t gene_start_idx) const {
    std::vector<std::string> out;
    for (size_t c = gene_start_idx; c < ncol(); ++c)
        out.push_back(col_names_[c]);
    return out;
}

} // namespace gsca
