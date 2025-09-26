#pragma once

#include "CSRMatrix.hpp"
#include "sparsityPattern.hpp"
#include "householderQR.hpp"
#include "mgsQR.hpp"
#include "eigenQRSolve.hpp"
#include "qr.hpp"

#include "launchThreads.hpp"

// @todo: add post filtration

template <typename T>
struct sam_internal_storage {
    std::vector<size_t> J;
    std::vector<T> rhs;
    std::vector<T> submatrix;
    std::vector<int> marker;
    QR<T> qr;

    sam_internal_storage(const size_t marker_size) : marker(marker_size, -1) {}
};

template <typename T, typename SparsityPatternType>
class SparseApproximateMap {
public:
    static void computeMap(const CSRMatrix<T>& targetMatrix, const CSRMatrix<T>& sourceMatrix, const SparsityPattern<T, SparsityPatternType>& sparsityPattern, CSRMatrix<T>& mappingMatrix) {
        // Get the computed sparsity pattern matrix
        const CSRMatrix<int> *pattern = sparsityPattern.getPattern();

        // Initialize the mapping matrix
        mappingMatrix.m_rows = targetMatrix.m_rows;
        mappingMatrix.m_cols = targetMatrix.m_cols;
        mappingMatrix.m_nnz = pattern->m_nnz;
        mappingMatrix.m_row_pointers = pattern->m_row_pointers;
        mappingMatrix.m_col_indices = pattern->m_col_indices;
        mappingMatrix.m_vals.resize(mappingMatrix.m_nnz);
        
        // Start SAM computation
        tbb::enumerable_thread_specific<sam_internal_storage<T>> local_storages(mappingMatrix.m_cols);
        tbb::parallel_for(tbb::blocked_range<size_t>(0, mappingMatrix.m_rows), [&](const tbb::blocked_range<size_t>& r) {
            sam_internal_storage<T>& ls = local_storages.local();
            for (size_t i = r.begin(); i < r.end(); ++i) {
                // Compute the submatrix indices for each row
                const size_t rowStart = pattern->m_row_pointers[i];
                const size_t rowEnd = pattern->m_row_pointers[i + 1];
                size_t iSize = rowEnd - rowStart;

                if (iSize == 0) continue;

                ls.J.clear();
                for (size_t j = rowStart; j < rowEnd; ++j) {
                    const size_t k = pattern->m_col_indices[j];

                    for (size_t k_start = pattern->m_row_pointers[k]; k_start < pattern->m_row_pointers[k + 1]; ++k_start) {
                        const size_t j_col = pattern->m_col_indices[k_start];
                        if (ls.marker[j_col] == -1) {
                            ls.marker[j_col] = 1;
                            ls.J.push_back(j_col);
                        }
                    }
                }
                std::sort(ls.J.begin(), ls.J.end());
                const size_t jSize = ls.J.size();

                // Use the marker as a map from sparse index to dense index (0 to |J|-1)
                for (size_t k = 0; k < jSize; ++k) {
                    ls.marker[ls.J[k]] = k;
                }

                // Compute the RHS vector from the target matrix
                ls.rhs.resize(jSize);
                std::fill(ls.rhs.begin(), ls.rhs.end(), T(0));
                for (size_t k = targetMatrix.m_row_pointers[i]; k < targetMatrix.m_row_pointers[i + 1]; ++k) {
                    const size_t col = targetMatrix.m_col_indices[k];
                    const int dense_idx = ls.marker[col];
                    if (dense_idx != -1) {
                        ls.rhs[dense_idx] = targetMatrix.m_vals[k];
                    }
                }

                // Compute the submatrix
                ls.submatrix.resize(jSize * iSize);
                std::fill(ls.submatrix.begin(), ls.submatrix.end(), T(0));
                for (size_t j = 0; j < iSize; ++j)
                {
                    const size_t source_row_idx = pattern->m_col_indices[rowStart + j];

                    for (size_t k = sourceMatrix.m_row_pointers[source_row_idx]; k < sourceMatrix.m_row_pointers[source_row_idx + 1]; ++k) {
                        const size_t col = sourceMatrix.m_col_indices[k];
                        const int dense_idx = ls.marker[col];
                        if (dense_idx != -1) {
                            ls.submatrix[dense_idx + jSize * j] = sourceMatrix.m_vals[k];
                        }
                    }
                }

                // Cleanup the marker for the next iteration
                for (size_t col : ls.J) {
                    ls.marker[col] = -1;
                }

                ls.qr.solve(jSize, iSize, ls.submatrix.data(), ls.rhs.data(), &mappingMatrix.m_vals[rowStart], storage_order::col_major);
            }
        });
    }

    static void post_filtration(const CSRMatrix<T>& map, CSRMatrix<T>& final_map, T threshold) {
        final_map.m_rows = map.m_rows;
        final_map.m_cols = map.m_cols;
        final_map.m_row_pointers.resize(final_map.m_rows + 1, 0);

        // count non-zero in the final map
        tbb::parallel_for(tbb::blocked_range<size_t>(0, final_map.m_rows), [&](const tbb::blocked_range<size_t> &r) {
            for (size_t i = r.begin(); i < r.end(); ++i) {
                int count = 0;
                for (size_t j = map.m_row_pointers[i], e = map.m_row_pointers[i + 1]; j < e; ++j) {
                    if (map.m_vals[j] >= threshold) {
                        count++;
                    }
                }

                final_map.m_row_pointers[i + 1] = count;
            }
        });

        final_map.m_nnz = final_map.scanRowSize();
        final_map.m_col_indices.resize(final_map.m_nnz, 0);
        final_map.m_vals.resize(final_map.m_nnz, T(0));

        tbb::parallel_for(tbb::blocked_range<size_t>(0, final_map.m_rows), [&](const tbb::blocked_range<size_t>& r) {
            for (size_t i = r.begin(); i < r.end(); ++i) {
                size_t current_pos = final_map.m_row_pointers[i];
                
                for (size_t j = map.m_row_pointers[i]; j < map.m_row_pointers[i + 1]; ++j) {
                    if (map.m_vals[j] >= threshold) {
                        final_map.m_col_indices[current_pos] = map.m_col_indices[j];
                        final_map.m_vals[current_pos] = map.m_vals[j];
                        current_pos++;
                    }
                }
            }
        });
    }
};