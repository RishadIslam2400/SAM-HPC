#pragma once

#include "CSRMatrix.hpp"
#include "launchThreads.hpp"

#include <atomic>

namespace sam {

void saad_pattern_extension(CSRMatrix<int> &prevPattern, CSRMatrix<int> &nextPattern) {
    nextPattern.m_rows = prevPattern.m_cols;
    nextPattern.m_cols = prevPattern.m_rows;
    nextPattern.m_row_pointers.resize(prevPattern.m_rows + 1);

    // Count the number of nonzero elements in each row of nextPattern
    tbb::enumerable_thread_specific<std::vector<int>> local_markers_1(nextPattern.m_cols, -1);
    tbb::parallel_for(tbb::blocked_range<size_t>(0, prevPattern.m_rows), [&](const tbb::blocked_range<size_t>& r) {
        std::vector<int>& marker = local_markers_1.local();
        for (size_t ia = r.begin(); ia < r.end(); ++ia) {
            int count = 0;
            const size_t rowBegA = prevPattern.m_row_pointers[ia];
            const size_t rowEndA = prevPattern.m_row_pointers[ia + 1];
            for (size_t ja = rowBegA; ja < rowEndA; ++ja) {
                const size_t colIdxA = prevPattern.m_col_indices[ja];
                const size_t rowBegB = prevPattern.m_row_pointers[colIdxA];
                const size_t rowEndB = prevPattern.m_row_pointers[colIdxA + 1];

                for (size_t jb = rowBegB; jb < rowEndB; ++jb) {
                    const size_t colIdxB = prevPattern.m_col_indices[jb];
                    if (marker[colIdxB] != static_cast<int>(ia)) {
                        marker[colIdxB] = static_cast<int>(ia);
                        ++count;
                    }
                }
            }
            nextPattern.m_row_pointers[ia + 1] = count;
        }
    });

    nextPattern.m_nnz = nextPattern.scanRowSize();
    nextPattern.m_col_indices.resize(nextPattern.m_nnz, 0);
    nextPattern.m_vals.resize(nextPattern.m_nnz, 1);

    // Compute the column indices
    local_markers_1.clear();
    // tbb::enumerable_thread_specific<std::vector<int>> local_markers_2(nextPattern.m_cols, -1);
    tbb::parallel_for(tbb::blocked_range<size_t>(0, nextPattern.m_rows), [&](const tbb::blocked_range<size_t>& r) {
        std::vector<int>& marker = local_markers_1.local();
        for (size_t ia = r.begin(); ia < r.end(); ++ia) {
            const size_t rowBeg = nextPattern.m_row_pointers[ia];
            size_t rowEnd = rowBeg;

            for (size_t ja = prevPattern.m_row_pointers[ia], ea = prevPattern.m_row_pointers[ia + 1]; ja < ea; ++ja) {
                size_t colIdxA = prevPattern.m_col_indices[ja];

                for (size_t jb = prevPattern.m_row_pointers[colIdxA], eb = prevPattern.m_row_pointers[colIdxA + 1]; jb < eb; ++jb) {
                    size_t colIdxB = prevPattern.m_col_indices[jb];
                    if (marker[colIdxB] < static_cast<int>(rowBeg)) {
                        marker[colIdxB] = static_cast<int>(rowEnd);
                        nextPattern.m_col_indices[rowEnd] = colIdxB;
                        ++rowEnd;
                    }
                }
            }

            for (size_t k = rowBeg; k < rowEnd; ++k) {
                marker[nextPattern.m_col_indices[k]] = -1;
            }

            std::sort(nextPattern.m_col_indices.begin() + rowBeg, nextPattern.m_col_indices.begin() + rowEnd);
        }
    });
}

template <bool needOut>
size_t* mergeRows(const size_t* row1, const size_t* row1_end, const size_t* row2, const size_t* row2_end, size_t* result) {
    while (row1 != row1_end && row2 != row2_end) {
        const size_t r1 = *row1;
        const size_t r2 = *row2;

        if (r1 < r2) {
            if constexpr (needOut) *result = r1;
            ++row1;
        } else if (r1 == r2) {
            if constexpr (needOut) *result = r1;
            ++row1;
            ++row2;
        } else {
            if constexpr (needOut) *result = r2;
            ++row2;
        }
        ++result;
    }

    if constexpr (needOut) {
        if (row1 < row1_end) {
            return std::copy(row1, row1_end, result);
        } else if (row2 < row2_end) {
            return std::copy(row2, row2_end, result);
        } else {
            return result;
        }
    } else {
        return result + (row1_end - row1) + (row2_end - row2);
    }
}

size_t prodRowWidth(const size_t *arow, const size_t *arow_end, const size_t *bptr,
                    const size_t *bcol, size_t *tmp_row1, size_t *tmp_row2, size_t *tmp_row3)
{
    const size_t nrows = arow_end - arow;

    // No rows merge, nothing to do
    if (nrows == 0)
        return 0;

    // Single row, just copy it to output
    if (nrows == 1)
        return bptr[*arow + 1] - bptr[*arow];

    // Two rows, merge them
    if (nrows == 2) {
        const size_t row1 = arow[0];
        const size_t row2 = arow[1];

        return mergeRows<false>(bcol + bptr[row1], bcol + bptr[row1 + 1],
                                bcol + bptr[row2], bcol + bptr[row2 + 1],
                                tmp_row1) - tmp_row1;
    }

    /**
     * Generic case (more than two rows).
     *
     * Merge rows by pairs, then merge the results together. When merging two rows, the result is always wider (or equal).
     * Merging by pairs allows to work with short rows as often as possible.
     */
    // merge first two rows
    size_t r1 = *arow++;
    size_t r2 = *arow++;

    size_t ncols1 = mergeRows<true>(bcol + bptr[r1], bcol + bptr[r1 + 1],
                                    bcol + bptr[r2], bcol + bptr[r2 + 1],
                                    tmp_row1) - tmp_row1;

    // Go by pairs
    while (arow + 1 < arow_end) {
        r1 = *arow++;
        r2 = *arow++;

        size_t ncols2 = mergeRows<true>(bcol + bptr[r1], bcol + bptr[r1 + 1],
                                        bcol + bptr[r2], bcol + bptr[r2 + 1],
                                        tmp_row2) - tmp_row2;

        if (arow == arow_end) {
            return mergeRows<false>(tmp_row1, tmp_row1 + ncols1,
                                    tmp_row2, tmp_row2 + ncols2,
                                    tmp_row3) - tmp_row3;
        } else {
            ncols1 = mergeRows<true>(tmp_row1, tmp_row1 + ncols1,
                                     tmp_row2, tmp_row2 + ncols2,
                                     tmp_row3) - tmp_row3;
            std::swap(tmp_row1, tmp_row3);
        }
    }

    // Merge the tail
    const size_t tail = *arow++;
    return mergeRows<false>(tmp_row1, tmp_row1 + ncols1,
                            bcol + bptr[tail], bcol + bptr[tail + 1],
                            tmp_row2) - tmp_row2;
}

void prodRow(const size_t* arow, const size_t* arow_end, const size_t* bptr, const size_t* bcol,
             size_t* out_row, size_t* tmp_row2, size_t* tmp_row3) {
    const size_t nrows = arow_end - arow;

    // No rows to merge, nothing to do
    if (nrows == 0)
        return;

    // Single row, just copy it to output
    if (nrows == 1) {
        const size_t idx = *arow;

        const size_t* browStart = bcol + bptr[idx];
        const size_t* browEnd = bcol + bptr[idx + 1];

        while (browStart != browEnd) {
            *out_row++ = *browStart++;
        }

        return;
    }

    // Two rows, merge them
    if (nrows == 2) {
        const size_t row_colind1 = arow[0];
        const size_t row_colind2 = arow[1];

        mergeRows<true>(bcol + bptr[row_colind1], bcol + bptr[row_colind1 + 1],
                        bcol + bptr[row_colind2], bcol + bptr[row_colind2 + 1], out_row);

        return;
    }

    /**
     * Generic case (more than two rows).
     *
     * Merge rows by pairs, then merge the results together. When merging two rows, the result is always wider (or equal).
     * Merging by pairs allows to work with short rows as often as possible.
     */
    // Merge first two rows
    size_t r1 = *arow++;
    size_t r2 = *arow++;

    size_t *tmp_row1 = out_row;

    size_t c_numcol1 = mergeRows<true>(bcol + bptr[r1], bcol + bptr[r1 + 1],
                                       bcol + bptr[r2], bcol + bptr[r2 + 1],
                                       tmp_row1) - tmp_row1;

    // Go by pairs
    while (arow + 1 < arow_end) {
        r1 = *arow++;
        r2 = *arow++;

        size_t c_numcol2 = mergeRows<true>(bcol + bptr[r1], bcol + bptr[r1 + 1],
                                           bcol + bptr[r2], bcol + bptr[r2 + 1],
                                           tmp_row2) - tmp_row2;

        c_numcol1 = mergeRows<true>(tmp_row1, tmp_row1 + c_numcol1,
                                    tmp_row2, tmp_row2 + c_numcol2,
                                    tmp_row3) - tmp_row3;

        std::swap(tmp_row3, tmp_row1);
    }

    // Merge the tail
    if (arow < arow_end) {
        r2 = *arow++;

        c_numcol1 = mergeRows<true>(tmp_row1, tmp_row1 + c_numcol1,
                                    bcol + bptr[r2], bcol + bptr[r2 + 1],
                                    tmp_row3) - tmp_row3;

        std::swap(tmp_row3, tmp_row1);
    }

    if (tmp_row1 != out_row) {
        std::copy(tmp_row1, tmp_row1 + c_numcol1, out_row);
    }
}

void rmerge_pattern_extension(CSRMatrix<int> &prevPattern, CSRMatrix<int> &nextPattern) {
    nextPattern.m_rows = prevPattern.m_rows;
    nextPattern.m_cols = prevPattern.m_cols;
    nextPattern.m_row_pointers.resize(nextPattern.m_rows + 1, 0);

    size_t maxRowWidth = tbb::parallel_reduce(
        tbb::blocked_range<size_t>(0, prevPattern.m_rows),
        size_t(0), // Identity for max operation
        [&](const tbb::blocked_range<size_t>& r, size_t localMax) -> size_t {
            for (size_t i = r.begin(); i < r.end(); ++i) {
                size_t rowWidth = 0;
                for (size_t j = prevPattern.m_row_pointers[i]; j < prevPattern.m_row_pointers[i + 1]; ++j) {
                    size_t colIdx = prevPattern.m_col_indices[j];
                    rowWidth += prevPattern.m_row_pointers[colIdx + 1] - prevPattern.m_row_pointers[colIdx];
                }
                localMax = std::max(localMax, rowWidth);
            }
            return localMax;
        },
        [](size_t a, size_t b) { return std::max(a, b); } // Join operation
    );

    

    // Temporary row of C for each thread.
    struct internal_storage {
        std::vector<size_t> tempCol;
        internal_storage(size_t maxWidth) {
            tempCol.resize(3 * maxWidth);
        }
    };
    tbb::enumerable_thread_specific<internal_storage> local_storages(maxRowWidth);

    tbb::parallel_for(tbb::blocked_range<size_t>(0, nextPattern.m_rows), 
        [&](const tbb::blocked_range<size_t>& r) {
            internal_storage& ls = local_storages.local();
            for (size_t i = r.begin(); i < r.end(); ++i) {
                const size_t rowStart = prevPattern.m_row_pointers[i];
                const size_t rowEnd = prevPattern.m_row_pointers[i + 1];
                nextPattern.m_row_pointers[i + 1] = prodRowWidth(
                    prevPattern.m_col_indices.data() + rowStart, prevPattern.m_col_indices.data() + rowEnd,
                    prevPattern.m_row_pointers.data(), prevPattern.m_col_indices.data(),
                    ls.tempCol.data(), ls.tempCol.data() + maxRowWidth, 
                    ls.tempCol.data() + 2 * maxRowWidth);
            }
        }
    );

    nextPattern.m_nnz = nextPattern.scanRowSize();
    nextPattern.m_col_indices.resize(nextPattern.m_nnz, 0);
    nextPattern.m_vals.resize(nextPattern.m_nnz, 1);

    tbb::parallel_for(tbb::blocked_range<size_t>(0, nextPattern.m_rows),
        [&](const tbb::blocked_range<size_t>& r) {
            internal_storage& ls = local_storages.local();
            for (size_t i = r.begin(); i < r.end(); ++i) {
                const size_t rowStart = prevPattern.m_row_pointers[i];
                const size_t rowEnd = prevPattern.m_row_pointers[i + 1];
                prodRow(
                    prevPattern.m_col_indices.data() + rowStart, prevPattern.m_col_indices.data() + rowEnd,
                    prevPattern.m_row_pointers.data(), prevPattern.m_col_indices.data(),
                    nextPattern.m_col_indices.data() + nextPattern.m_row_pointers[i],
                    ls.tempCol.data(), ls.tempCol.data() + maxRowWidth);
            }
        }
    );
}

void extend_pattern(CSRMatrix<int> &S, const int level) {

    std::vector<CSRMatrix<int>> finalPatterns(level - 1);

    for (int l = 0; l < level - 1; ++l) {
        if (l == 0) {
            if (num_threads < 16) {
                saad_pattern_extension(S, finalPatterns[l]);
            } else {
                rmerge_pattern_extension(S, finalPatterns[l]);
            }
        } else if (l >= 1) {
            if (num_threads < 16) {
                saad_pattern_extension(finalPatterns[l - 1], finalPatterns[l]);
            } else {
                rmerge_pattern_extension(finalPatterns[l - 1], finalPatterns[l]);
            }
        }
    }

    swap(S, finalPatterns.back());
}

} // namespace sam