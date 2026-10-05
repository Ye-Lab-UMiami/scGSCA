// ============================================================================
// main.cpp
//
// Batch runner for MetaGSCA analysis on KEGG pathway files.
// This replaces examples_all_kegg.R.
//
// Usage:
//   scgsca_run --kegg-dir <dir> --out-dir <dir> [options]
//
// Each TSV file in kegg-dir is one KEGG pathway.
// DC and CT tests can be enabled or disabled independently.
// ============================================================================
#include "MetaGSCA.h"
#include <algorithm>
#include <chrono>
#include <cctype>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <set>
#include <sstream>

namespace fs = std::filesystem;

struct CTPair {
    std::string ct1;
    std::string ct2;
    std::string suffix;
};

// ---------------------------------------------------------------------------
// Parse command-line arguments
// ---------------------------------------------------------------------------
struct CLIArgs {
    std::string kegg_dir;
    std::string out_dir;
    int nperm       = 100;
    uint64_t seed   = 416627;
    int cores       = 8;
    char ref_letter = 'A';
    int gene_start  = 3;           // 0-based column index where genes begin
    bool skip_dc    = false;
    bool skip_ct    = false;
    bool save_ct_perm_corr = false;
    bool paired = false;
    std::vector<int> file_indices; // 1-based; empty = all
    std::vector<CTPair> ct_pairs;
};

static std::string normalize_token(const std::string& s) {
    std::string out;
    out.reserve(s.size());
    for (unsigned char ch : s) {
        if (std::isalnum(ch))
            out.push_back(static_cast<char>(std::tolower(ch)));
    }
    return out;
}

static std::string sanitize_suffix(const std::string& s) {
    std::string out;
    bool prev_us = false;
    for (unsigned char ch : s) {
        if (std::isalnum(ch)) {
            out.push_back(static_cast<char>(std::tolower(ch)));
            prev_us = false;
        } else if (!prev_us && !out.empty()) {
            out.push_back('_');
            prev_us = true;
        }
    }
    while (!out.empty() && out.back() == '_')
        out.pop_back();
    return out.empty() ? "pair" : out;
}

static CTPair parse_ct_pair(const std::string& spec) {
    auto pos = spec.find(',');
    if (pos == std::string::npos)
        throw std::runtime_error(
            "Invalid --ct-pair value '" + spec + "'. Expected ct1,ct2.");

    CTPair pair;
    pair.ct1 = spec.substr(0, pos);
    pair.ct2 = spec.substr(pos + 1);
    if (pair.ct1.empty() || pair.ct2.empty())
        throw std::runtime_error(
            "Invalid --ct-pair value '" + spec + "'. Empty cell type name.");
    pair.suffix = sanitize_suffix(pair.ct1) + "_vs_" + sanitize_suffix(pair.ct2);
    return pair;
}

static std::string find_matching_col(const gsca::DataFrame& dt,
                                     const std::vector<std::string>& aliases) {
    std::set<std::string> alias_norm;
    for (const auto& alias : aliases)
        alias_norm.insert(normalize_token(alias));

    for (const auto& name : dt.col_names()) {
        if (alias_norm.count(normalize_token(name)))
            return name;
    }
    return "";
}

static void normalize_input_columns(gsca::DataFrame& dt) {
    const struct AliasMap {
        std::string canonical;
        std::vector<std::string> aliases;
    } maps[] = {
        {"Cell_ID", {"Cell_ID", "cell_id", "CellID", "cellid"}},
        {"Celltype", {"Celltype", "cell_type", "celltype"}},
        {"Samples", {"Samples", "samples", "sample"}},
    };

    for (const auto& mapping : maps) {
        if (dt.has_col(mapping.canonical))
            continue;
        auto found = find_matching_col(dt, mapping.aliases);
        if (!found.empty() && found != mapping.canonical)
            dt.rename_col(found, mapping.canonical);
    }
}

static bool pair_is_present_in_both_groups(const gsca::DataFrame& dt,
                                           const std::vector<int>& grp,
                                           const CTPair& pair) {
    const auto& ct = dt.str_col("Celltype");
    bool g1_ct1 = false, g1_ct2 = false, g2_ct1 = false, g2_ct2 = false;
    for (size_t i = 0; i < dt.nrow(); ++i) {
        if (grp[i] == 1 && ct[i] == pair.ct1) g1_ct1 = true;
        if (grp[i] == 1 && ct[i] == pair.ct2) g1_ct2 = true;
        if (grp[i] == 2 && ct[i] == pair.ct1) g2_ct1 = true;
        if (grp[i] == 2 && ct[i] == pair.ct2) g2_ct2 = true;
    }
    return g1_ct1 && g1_ct2 && g2_ct1 && g2_ct2;
}

static bool run_ct_pair(const std::string& tag,
                        const CTPair& pair,
                        const gsca::GSCAConfig& base_cfg,
                        gsca::DataFrame& dt,
                        const std::vector<int>& grp,
                        char ref_letter,
                        const std::vector<std::string>& genelist) {
    std::vector<bool> mask(dt.nrow(), false);
    const auto& ct = dt.str_col("Celltype");
    size_t kept = 0;
    for (size_t i = 0; i < dt.nrow(); ++i) {
        if (ct[i] == pair.ct1 || ct[i] == pair.ct2) {
            mask[i] = true;
            ++kept;
        }
    }

    if (kept == 0) {
        std::cout << "Skipping CT pair " << pair.ct1 << " vs " << pair.ct2
                  << " for " << tag << ": neither cell type is present.\n";
        return false;
    }

    auto dt_pair = dt.row_subset(mask);
    std::vector<int> grp_pair;
    grp_pair.reserve(kept);
    for (size_t i = 0; i < dt.nrow(); ++i)
        if (mask[i]) grp_pair.push_back(grp[i]);

    if (!pair_is_present_in_both_groups(dt_pair, grp_pair, pair)) {
        std::cout << "Skipping CT pair " << pair.ct1 << " vs " << pair.ct2
                  << " for " << tag
                  << ": one or both disease groups do not contain both cell types.\n";
        return false;
    }

    auto cfg = base_cfg;
    cfg.projectname = tag + ".ct." + pair.suffix;
    cfg.ct_ref_level = pair.ct1;
    gsca::GSAR_qr_ct(cfg, dt_pair, grp_pair, static_cast<int>(ref_letter), genelist);
    return true;
}

static CLIArgs parse_args(int argc, char* argv[]) {
    CLIArgs args;
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--kegg-dir" && i + 1 < argc) args.kegg_dir = argv[++i];
        else if (arg == "--out-dir" && i + 1 < argc) args.out_dir = argv[++i];
        else if (arg == "--nperm" && i + 1 < argc) args.nperm = std::stoi(argv[++i]);
        else if (arg == "--seed" && i + 1 < argc) args.seed = std::stoull(argv[++i]);
        else if (arg == "--cores" && i + 1 < argc) args.cores = std::stoi(argv[++i]);
        else if (arg == "--ref-letter" && i + 1 < argc) args.ref_letter = argv[++i][0];
        else if (arg == "--gene-start" && i + 1 < argc) args.gene_start = std::stoi(argv[++i]);
        else if (arg == "--skip-dc") args.skip_dc = true;
        else if (arg == "--skip-ct") args.skip_ct = true;
        else if (arg == "--save-ct-perm-corr") args.save_ct_perm_corr = true;
        else if (arg == "--paired") args.paired = true;
        else if (arg == "--unpaired") args.paired = false;
        else if (arg == "--ct-tail" || arg.rfind("--ct-tail=", 0) == 0) {
            throw std::invalid_argument("--ct-tail is not supported. CT outputs both absolute-delta (two-sided) and signed right-tail P values.");
        }
        else if (arg == "--ct-pair" && i + 1 < argc) {
            args.ct_pairs.push_back(parse_ct_pair(argv[++i]));
        }
        else if (arg == "--files" && i + 1 < argc) {
            std::string list = argv[++i];
            std::istringstream ss(list);
            std::string tok;
            while (std::getline(ss, tok, ','))
                args.file_indices.push_back(std::stoi(tok));
        }
        else if (arg == "--help" || arg == "-h") {
            std::cout
                << "Usage: scgsca_run --kegg-dir <dir> --out-dir <dir> [options]\n"
                << "  --kegg-dir <path>    KEGG .tsv file directory\n"
                << "  --out-dir <path>     Output root directory\n"
                << "  --nperm <int>        Number of permutations (default: 100)\n"
                << "  --seed <int>         Random seed (default: 416627)\n"
                << "  --cores <int>        CPU cores (default: 8)\n"
                << "  --ref-letter <char>  Reference group letter (default: A)\n"
                << "  --gene-start <int>   0-based gene column index (default: 3)\n"
                << "  --files <idx,...>    1-based file indices to run (default: all)\n"
                << "  --skip-dc            Skip the sample-phenotype test.\n"
                << "  --paired             Use the paired subject-level permutation procedure.\n"
                << "  --unpaired           Use unpaired subject-label permutations (default).\n"
                << "  --save-ct-perm-corr  Save CT permutation correlation matrices as TSVs.\n"
                << "  --ct-pair a,b        Run the cell-type interaction test for this pair.\n"
                << "                       Repeat to run multiple CT comparisons.\n"
                << "  --skip-ct            Skip the cell-type interaction test.\n";
            std::cout << "  CT test: two-sided absolute-delta P value: Permutation.Pvalue.AbsDelta.\n"
                      << "           Additional signed output: fixed right tail (Permutation.Pvalue.Delta).\n"
                      << "           Delta = D(non-reference group) - D(reference group).\n"
                      << "           Select the reference with --ref-letter.\n";
            std::exit(0);
        }
    }

    if (args.ct_pairs.empty() && !args.skip_ct) {
        args.ct_pairs.push_back({"Astrocyte", "Endothelial_cells", "AvE"});
        args.ct_pairs.push_back({"Astrocyte", "Macrophage", "AvM"});
    }

    return args;
}

// ---------------------------------------------------------------------------
// Run one pathway: DC + optional CT tests
// ---------------------------------------------------------------------------
static bool run_one_test(const std::string& tag,
                         gsca::DataFrame& dt,
                         char ref_letter,
                         const std::string& resdir,
                         int gene_start,
                         int nperm,
                         uint64_t seed,
                         int cores,
                         const std::vector<CTPair>& ct_pairs,
                         bool save_ct_perm_corr,
                         bool skip_dc,
                         bool skip_ct,
                         bool paired)
{
    normalize_input_columns(dt);

    // Validate required columns
    if (!dt.has_col("Samples"))
        throw std::runtime_error(tag + ": column 'Samples' not found.");
    if (!dt.has_col("Celltype"))
        throw std::runtime_error(tag + ": column 'Celltype' not found.");

    // Determine group from last character of Samples column
    const auto& samples = dt.str_col("Samples");
    std::vector<int> grp(samples.size());
    std::set<char> suffixes;
    for (size_t i = 0; i < samples.size(); ++i) {
        if (samples[i].empty())
            throw std::runtime_error(tag + ": empty value found in Samples column.");
        char last = samples[i].back();
        suffixes.insert(last);
    }

    if (suffixes.size() != 2)
        throw std::runtime_error(
            tag + ": Samples suffix has " + std::to_string(suffixes.size()) +
            " levels. Expected exactly 2.");

    if (suffixes.find(ref_letter) == suffixes.end())
        throw std::runtime_error(
            tag + ": ref_letter '" + std::string(1, ref_letter) +
            "' not found in Samples suffix levels.");

    for (size_t i = 0; i < samples.size(); ++i)
        grp[i] = (samples[i].back() == ref_letter) ? 1 : 2;

    std::cout << "\n============================================================\n"
              << "Running pathway: " << tag
              << " (group.ref = " << ref_letter << ")\n"
              << "Output folder: " << resdir << "\n"
              << "nperm = " << nperm << "  seed = " << seed
              << "  cores = " << cores
              << "  paired = " << (paired ? "TRUE" : "FALSE") << "\n"
              << "============================================================\n\n";

    gsca::GSCAConfig cfg;
    cfg.gene_start = static_cast<size_t>(gene_start);
    cfg.ct_var     = "Celltype";
    cfg.subj_var   = "Samples";
    cfg.nperm      = nperm;
    cfg.seed       = seed;
    cfg.num_cores  = cores;
    cfg.check_sd   = true;
    cfg.min_sd     = 0.001;
    cfg.filter_count = 0.0;
    cfg.na_counts  = true;
    cfg.cor_method = "pearson";
    cfg.max_skip   = 50;
    cfg.n_sim      = 100;
    cfg.diagnostic_plots = false;  // not implemented in C++
    cfg.save_ct_perm_corr = save_ct_perm_corr;
    cfg.paired      = paired;
    cfg.resdir      = resdir;

    std::vector<std::string> genelist = {"NA"};  // use all genes

    if (skip_dc) {
        std::cout << "Skipping DC test because --skip-dc was requested.\n";
    } else {
        cfg.projectname = tag + ".dc";
        gsca::GSAR_qr_dc(cfg, dt, grp, static_cast<int>(ref_letter), genelist);
    }

    if (skip_ct) {
        std::cout << "Skipping all CT tests because --skip-ct was requested.\n";
        return true;
    }

    size_t n_ct_run = 0;
    char other_letter = *std::find_if(suffixes.begin(), suffixes.end(),
                                    [ref_letter](char c) { return c != ref_letter; });
    std::cout << "CT two-sided output: Permutation.Pvalue.AbsDelta.\n";
    std::cout << "CT tail: right; Delta = D(" << other_letter
              << ") - D(" << ref_letter << "). Applies only to the additional signed-delta output.\n";
    for (const auto& pair : ct_pairs) {
        if (run_ct_pair(tag, pair, cfg, dt, grp, ref_letter, genelist))
            ++n_ct_run;
    }

    if (n_ct_run == 0)
        std::cout << "No CT comparisons were run for " << tag << ".\n";

    return true;
}

// =========================================================================
// main
// =========================================================================
int main(int argc, char* argv[]) {
    CLIArgs args;
    try { args = parse_args(argc, argv); }
    catch (const std::exception& e) {
        std::cerr << "Argument error: " << e.what() << "\n";
        return 1;
    }

    if (args.kegg_dir.empty() || args.out_dir.empty()) {
        std::cerr << "Error: --kegg-dir and --out-dir are required.\n"
                  << "Use --help for usage.\n";
        return 1;
    }
    if (args.skip_dc && args.skip_ct) {
        std::cerr << "Error: both --skip-dc and --skip-ct were requested; nothing to run.\n";
        return 1;
    }

    // Collect .tsv files
    std::vector<std::string> all_files;
    for (auto& entry : fs::directory_iterator(args.kegg_dir)) {
        if (entry.path().extension() == ".tsv")
            all_files.push_back(entry.path().filename().string());
    }
    std::sort(all_files.begin(), all_files.end());

    std::cout << "KEGG directory: " << args.kegg_dir << "\n"
              << "Found " << all_files.size() << " TSV files:\n";
    for (size_t i = 0; i < all_files.size(); ++i)
        std::cout << "  [" << (i + 1) << "] " << all_files[i] << "\n";

    if (all_files.empty()) {
        std::cerr << "No .tsv files found in " << args.kegg_dir << ".\n";
        return 1;
    }

    if (args.skip_dc) {
        std::cout << "DC mode: disabled (--skip-dc)\n";
    } else {
        std::cout << "DC mode: enabled\n";
    }

    if (args.skip_ct) {
        std::cout << "CT mode: disabled (--skip-ct)\n";
    } else {
        std::cout << "CT pairs to try:\n";
        for (const auto& pair : args.ct_pairs) {
            std::cout << "  - " << pair.ct1 << " vs " << pair.ct2
                      << " -> suffix " << pair.suffix << "\n";
        }
    }

    // Select files to run
    std::vector<std::string> files_to_run;
    if (args.file_indices.empty()) {
        files_to_run = all_files;
    } else {
        for (int idx : args.file_indices) {
            if (idx >= 1 && idx <= static_cast<int>(all_files.size()))
                files_to_run.push_back(all_files[idx - 1]);
        }
    }
    std::cout << "\nFiles selected to run: " << files_to_run.size() << "\n\n";

    // Create output directory
    fs::create_directories(args.out_dir);

    // Batch loop
    struct Summary {
        std::string pathway;
        bool success;
        double time_sec;
    };
    std::vector<Summary> results;

    for (auto& fname : files_to_run) {
        std::string fpath = args.kegg_dir + "/" + fname;
        std::string pathway = fname.substr(0, fname.size() - 4); // strip .tsv
        std::string out_path = args.out_dir + "/" + pathway;
        fs::create_directories(out_path);

        std::cout << "Reading file: " << fpath << "\n";

        auto t_start = std::chrono::steady_clock::now();
        bool ok = false;

        try {
            auto dt = gsca::read_tsv(fpath);
            normalize_input_columns(dt);

            std::cout << "Rows: " << dt.nrow() << "\n"
                      << "Columns: " << dt.ncol() << "\n";

            std::cout << "Column names: ";
            for (size_t i = 0; i < dt.ncol(); ++i) {
                if (i > 0) std::cout << ", ";
                std::cout << "[" << dt.col_name(i) << "]";
            }
            std::cout << "\n";
            std::cout << "Has 'Samples': " << (dt.has_col("Samples") ? "YES" : "NO") << "\n";
            std::cout << "Has 'Celltype': " << (dt.has_col("Celltype") ? "YES" : "NO") << "\n\n";

            ok = run_one_test(pathway, dt, args.ref_letter, out_path,
                              args.gene_start, args.nperm, args.seed,
                              args.cores, args.ct_pairs, args.save_ct_perm_corr,
                              args.skip_dc, args.skip_ct, args.paired);
        } catch (std::exception& e) {
            std::cerr << "\nERROR for pathway: " << pathway << "\n"
                      << e.what() << "\n";
            ok = false;
        }

        auto t_end = std::chrono::steady_clock::now();
        double elapsed = std::chrono::duration<double>(t_end - t_start).count();

        std::cout << "\nFinished: " << pathway
                  << " | success = " << (ok ? "TRUE" : "FALSE")
                  << " | time = " << elapsed << " sec\n";

        results.push_back({pathway, ok, elapsed});
    }

    std::cout << "\n\n========== BATCH SUMMARY ==========\n";
    int n_success = 0, n_fail = 0;
    for (auto& r : results) {
        std::cout << "  " << r.pathway
                  << " | " << (r.success ? "OK" : "FAIL")
                  << " | " << r.time_sec << " sec\n";
        if (r.success) ++n_success; else ++n_fail;
    }
    std::cout << "Total: " << results.size()
              << " | Succeeded: " << n_success
              << " | Failed: " << n_fail << "\n";

    {
        std::ofstream ofs(args.out_dir + "/batch_summary.csv");
        ofs << "pathway,success,time_sec\n";
        for (auto& r : results)
            ofs << r.pathway << "," << (r.success ? "TRUE" : "FALSE")
                << "," << r.time_sec << "\n";
    }
    std::cout << "Summary saved to: " << args.out_dir << "/batch_summary.csv\n";

    return 0;
}
