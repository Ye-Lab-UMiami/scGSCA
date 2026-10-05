// ============================================================================
// gsnca_check.cpp
//
// Input validation for GSNCA analysis.
// Ports: gsnca_check_dc.R, gsnca_check_ct.R
// ============================================================================
#include "MetaGSCA.h"
#include <algorithm>
#include <set>

namespace gsca {

// ---------------------------------------------------------------------------
// gsnca_check_dc
//
// Validates inputs for DC (disease condition) mode:
//   - projectname must be non-empty
//   - group must have same length as object rows
//   - group values must be 1 or 2
//   - at least 3 samples per group
//   - group_ref must be specified (non-zero)
// ---------------------------------------------------------------------------
void gsnca_check_dc(const std::string& projectname,
                    const DataFrame& object,
                    const std::vector<int>& group,
                    int group_ref,
                    const std::string& ct_var)
{
    if (projectname.empty())
        throw std::runtime_error("You must provide a project name.");

    size_t nv = object.nrow();

    if (group.size() != nv)
        throw std::runtime_error(
            "length of 'group' must equal the number of rows in 'object'");

    int count1 = 0, count2 = 0;
    for (int g : group) {
        if (g != 1 && g != 2)
            throw std::runtime_error(
                "all members in 'group' must have values 1 or 2");
        if (g == 1) ++count1;
        if (g == 2) ++count2;
    }

    if (count1 < 3 || count2 < 3)
        throw std::runtime_error(
            "there are less than 3 samples in at least one group");

    if (group_ref == 0)
        throw std::runtime_error("Reference group must be specified.");
}

// ---------------------------------------------------------------------------
// gsnca_check_ct
//
// Validates inputs for CT (cell type) mode.
// Same as DC checks plus: at most 2 unique cell types.
// ---------------------------------------------------------------------------
void gsnca_check_ct(const std::string& projectname,
                    const DataFrame& object,
                    const std::vector<int>& group,
                    int group_ref,
                    const std::string& ct_var)
{
    // Run all DC checks first
    gsnca_check_dc(projectname, object, group, group_ref, ct_var);

    // Check that there are at most 2 unique cell types
    const auto& ct_col = object.str_col(ct_var);
    std::set<std::string> unique_ct(ct_col.begin(), ct_col.end());

    if (unique_ct.size() > 2)
        throw std::runtime_error(
            "We can only test 2 cell types. "
            "Subset data to the two cell types of interest for testing.");
}

} // namespace gsca
