// ============================================================================
// tsv_io.h
//
// Read / write tab-separated files into/from DataFrame.
// Mirrors R's data.table::fread and write.table with sep='\t'.
// ============================================================================
#ifndef METAGSCA_TSV_IO_H
#define METAGSCA_TSV_IO_H

#include "DataFrame.h"
#include <string>

namespace gsca {

// Read a TSV file.  Columns that are entirely parseable as double are stored
// as numeric; otherwise they are stored as string columns.
DataFrame read_tsv(const std::string& path);

// Write a DataFrame to a TSV file.
void write_tsv(const DataFrame& df, const std::string& path);

} // namespace gsca

#endif // METAGSCA_TSV_IO_H
