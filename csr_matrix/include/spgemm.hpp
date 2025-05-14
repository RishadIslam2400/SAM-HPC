#pragma once

#include "CSRMatrix.hpp"
#include "launchThreads.hpp"

#include <atomic>

template <typename T>
void spgemm_saad(const CSRMatrix<T> &A, const CSRMatrix<T> &B, CSRMatrix<T> &C)
{
    C.row_num = A.row_num;
    C.col_num = B.col_num;
    C.row_pointers = new std::vector<size_t>(C.row_num + 1, 0);

    // Count the number of non-zero elements in each row of C
    if constexpr (SEQUENTIAL) {
        std::vector<int> marker(C.col_num, -1);
        for (size_t ia = 0; ia < C.row_num; ++ia) {
            int count = 0;
            const size_t rowBegA = (*(A.row_pointers))[ia];
            const size_t rowEndA = (*(A.row_pointers))[ia + 1];
            for (size_t ja = rowBegA; ja < rowEndA; ++ja) {
                const size_t colIdxA = (*(A.col_indices))[ja];

                const size_t rowBegB = (*(B.row_pointers))[colIdxA];
                const size_t rowEndB = (*(B.row_pointers))[colIdxA + 1];
                for (size_t jb = rowBegB; jb < rowEndB; ++jb) {
                    const size_t colIdxB = (*(B.col_indices))[jb];
                    if (marker[colIdxB] != static_cast<int>(ia)) {
                        marker[colIdxB] = static_cast<int>(ia);
                        ++count;
                    }
                }
            }
            (*(C.row_pointers))[ia + 1] = count;
        }
    } else {
        auto f = [&](size_t start, size_t end) {
            std::vector<int> marker(C.col_num, -1);
            for (size_t ia = start; ia < end; ++ia) {
                int count = 0;
                const size_t rowBegA = (*(A.row_pointers))[ia];
                const size_t rowEndA = (*(A.row_pointers))[ia + 1];
                for (size_t ja = rowBegA; ja < rowEndA; ++ja) {
                    const size_t colIdxA = (*(A.col_indices))[ja];

                    const size_t rowBegB = (*(B.row_pointers))[colIdxA];
                    const size_t rowEndB = (*(B.row_pointers))[colIdxA + 1];
                    for (size_t jb = rowBegB; jb < rowEndB; ++jb) {
                        const size_t colIdxB = (*(B.col_indices))[jb];
                        if (marker[colIdxB] != static_cast<int>(ia)) {
                            marker[colIdxB] = static_cast<int>(ia);
                            ++count;
                        }
                    }
                }
                (*(C.row_pointers))[ia + 1] = count;
            }
        };

        launchThreads(C.row_num, f);
    }

    C.nnz = C.scanRowSize();
    C.col_indices = new std::vector<size_t>(C.nnz, 0);
    C.vals = new std::vector<T>(C.nnz, T());

    // Todo: remove span, use pointers
    // Compute the column indices and values of C
    if constexpr (SEQUENTIAL) {
        std::vector<int> marker(C.col_num, -1);
        for (size_t ia = 0; ia < C.row_num; ++ia) {
            const size_t rowBeg = (*(C.row_pointers))[ia];
            size_t rowEnd = rowBeg;
            for (size_t ja = (*(A.row_pointers))[ia], ea = (*(A.row_pointers))[ia + 1]; ja < ea; ++ja) {
                size_t colIdxA = (*(A.col_indices))[ja];
                T valA = (*(A.vals))[ja];

                for (size_t jb = (*(B.row_pointers))[colIdxA], eb = (*(B.row_pointers))[colIdxA + 1]; jb < eb; ++jb) {
                    size_t colIdxB = (*(B.col_indices))[jb];
                    T valB = (*(B.vals))[jb];
                    if (marker[colIdxB] < static_cast<int>(rowBeg)) {
                        marker[colIdxB] = static_cast<int>(rowEnd);
                        (*(C.col_indices))[rowEnd] = colIdxB;
                        (*(C.vals))[rowEnd] = valA * valB;
                        ++rowEnd;
                    }
                    else {
                        (*(C.vals))[marker[colIdxB]] += valA * valB;
                    }
                }
            }

            std::span<size_t> col_indices_span(C.col_indices->begin() + rowBeg, rowEnd - rowBeg);
            std::span<T> vals_span(C.vals->begin() + rowBeg, rowEnd - rowBeg);
            CSRMatrix<T>::sortRow(col_indices_span, vals_span, rowEnd - rowBeg);
        }
    } else {
        auto f = [&](size_t start, size_t end) {
            std::vector<int> marker(C.col_num, -1);
            for (size_t ia = start; ia < end; ++ia) {
                const size_t rowBeg = (*(C.row_pointers))[ia];
                size_t rowEnd = rowBeg;
                for (size_t ja = (*(A.row_pointers))[ia], ea = (*(A.row_pointers))[ia + 1]; ja < ea; ++ja) {
                    size_t colIdxA = (*(A.col_indices))[ja];
                    T valA = (*(A.vals))[ja];

                    for (size_t jb = (*(B.row_pointers))[colIdxA], eb = (*(B.row_pointers))[colIdxA + 1]; jb < eb; ++jb) {
                        size_t colIdxB = (*(B.col_indices))[jb];
                        T valB = (*(B.vals))[jb];
                        if (marker[colIdxB] < static_cast<int>(rowBeg)) {
                            marker[colIdxB] = static_cast<int>(rowEnd);
                            (*(C.col_indices))[rowEnd] = colIdxB;
                            (*(C.vals))[rowEnd] = valA * valB;
                            ++rowEnd;
                        }
                        else {
                            (*(C.vals))[marker[colIdxB]] += valA * valB;
                        }
                    }
                }

                std::span<size_t> col_indices_span(C.col_indices->begin() + rowBeg, rowEnd - rowBeg);
                std::span<T> vals_span(C.vals->begin() + rowBeg, rowEnd - rowBeg);
                CSRMatrix<T>::sortRow(col_indices_span, vals_span, rowEnd - rowBeg);
            }
        };

        launchThreads(C.row_num, f);
    }
}

template <bool needOut>
size_t *mergeRows(std::span<const size_t> row1, std::span<const size_t> row2, size_t *result)
{
    auto row1_it = row1.begin();
    auto row2_it = row2.begin();

    while (row1_it != row1.end() && row2_it != row2.end())
    {
        if (*row1_it < *row2_it)
        {
            if constexpr (needOut)
            {
                *result = *row1_it;
            }
            ++row1_it;
        }
        else if (*row1_it == *row2_it)
        {
            if constexpr (needOut)
            {
                *result = *row1_it;
            }
            ++row1_it;
            ++row2_it;
        }
        else
        {
            if constexpr (needOut)
            {
                *result = *row2_it;
            }
            ++row2_it;
        }
        ++result;
    }

    if constexpr (needOut)
    {
        if (row1_it < row1.end())
        {
            return std::copy(row1_it, row1.end(), result);
        }
        else if (row2_it < row2.end())
        {
            return std::copy(row2_it, row2.end(), result);
        }
        else
        {
            return result;
        }
    }
    else
    {
        return result + (row1.end() - row1_it) + (row2.end() - row2_it);
    }
}

template <typename T>
size_t *mergeRows(T alpha1, std::span<const size_t> row_colind1, std::span<const T> row_vals1,
                  T alpha2, std::span<const size_t> row_colind2, std::span<const T> row_vals2,
                  size_t *crow_colind, T *crow_vals)
{
    auto row_colind1_it = row_colind1.begin();
    auto row_colind2_it = row_colind2.begin();
    auto row_vals1_it = row_vals1.begin();
    auto row_vals2_it = row_vals2.begin();

    while (row_colind1_it != row_colind1.end() && row_colind2_it != row_colind2.end())
    {
        const size_t col1 = *row_colind1_it;
        const size_t col2 = *row_colind2_it;

        if (col1 < col2)
        {
            ++row_colind1_it;
            *crow_colind = col1;
            *crow_vals = alpha1 * (*row_vals1_it++);
        }
        else if (col1 == col2)
        {
            ++row_colind1_it;
            ++row_colind2_it;

            *crow_colind = col1;
            *crow_vals = alpha1 * (*row_vals1_it++) + alpha2 * (*row_vals2_it++);
        }
        else
        {
            ++row_colind2_it;
            *crow_colind = col2;
            *crow_vals = alpha2 * (*row_vals2_it++);
        }

        ++crow_colind;
        ++crow_vals;
    }

    while (row_colind1_it < row_colind1.end())
    {
        *crow_colind++ = *row_colind1_it++;
        *crow_vals++ = alpha1 * (*row_vals1_it++);
    }

    while (row_colind2_it < row_colind2.end())
    {
        *crow_colind++ = *row_colind2_it++;
        *crow_vals++ = alpha2 * (*row_vals2_it++);
    }

    return crow_colind;
}

size_t prodRowWidth(std::span<const size_t> arow, const std::vector<size_t> &browptr, const std::vector<size_t> &bcolind,
                    size_t *temp_col1, size_t *temp_col2, size_t *temp_col3)
{
    const size_t nrows = arow.size();

    // No rows merge, nothing to do
    if (nrows == 0)
        return 0;

    // Single row, just copy it to output
    if (nrows == 1)
    {
        return browptr[arow[0] + 1] - browptr[arow[0]];
    }

    // Two rows, merge them
    if (nrows == 2)
    {
        const size_t row1 = arow[0];
        const size_t row2 = arow[1];

        std::span<const size_t> brow1(bcolind.begin() + browptr[row1], browptr[row1 + 1] - browptr[row1]);
        std::span<const size_t> brow2(bcolind.begin() + browptr[row2], browptr[row2 + 1] - browptr[row2]);

        return mergeRows<false>(brow1, brow2, temp_col1) - temp_col1;
    }

    /**
     * Generic case (more than two rows).
     *
     * Merge rows by pairs, then merge the results together. When merging two rows, the result is always wider (or equal).
     * Merging by pairs allows to work with short rows as often as possible.
     */
    auto arow_it = arow.begin();
    // merge first two rows
    const auto r1 = *arow_it++;
    std::span<const size_t> brow1(bcolind.begin() + browptr[r1], browptr[r1 + 1] - browptr[r1]);
    const auto r2 = *arow_it++;
    std::span<const size_t> brow2(bcolind.begin() + browptr[r2], browptr[r2 + 1] - browptr[r2]);

    size_t ncols1 = mergeRows<true>(brow1, brow2, temp_col1) - temp_col1;

    // Go by pairs
    while (arow_it + 1 < arow.end())
    {
        const auto a1 = *arow_it++;
        const auto a2 = *arow_it++;
        std::span<const size_t> browfirst(bcolind.begin() + browptr[a1], browptr[a1 + 1] - browptr[a1]);
        std::span<const size_t> browsecond(bcolind.begin() + browptr[a2], browptr[a2 + 1] - browptr[a2]);

        size_t ncols2 = mergeRows<true>(browfirst, browsecond, temp_col2) - temp_col2;

        std::span<const size_t> temp_col1_span(temp_col1, ncols1);
        std::span<const size_t> temp_col2_span(temp_col2, ncols2);

        if (arow_it == arow.end())
        {
            return mergeRows<false>(temp_col1_span, temp_col2_span, temp_col3) - temp_col3;
        }
        else
        {
            ncols1 = mergeRows<true>(temp_col1_span, temp_col2_span, temp_col3) - temp_col3;
            std::swap(temp_col1, temp_col3);
        }
    }

    // Merge the tail
    const auto tail = *arow_it++;
    std::span<const size_t> browtail(bcolind.begin() + browptr[tail], browptr[tail + 1] - browptr[tail]);
    std::span<const size_t> temp_col1_span(temp_col1, ncols1);
    return mergeRows<false>(temp_col1_span, browtail, temp_col2) - temp_col2;
}

template <typename T>
void prodRow(std::span<const size_t> arow_colind, std::span<const T> arow_vals,
             const std::vector<size_t> &browptr, const std::vector<size_t> &bcolind, const std::vector<T> &bvals,
             size_t *crow_colind, T *crow_vals,
             size_t *temp_col2, T *temp_val2,
             size_t *temp_col3, T *temp_val3)
{
    const size_t nrows = arow_colind.size();

    // No rows to merge, nothing to do
    if (nrows == 0)
        return;

    // Single row, just copy it to output
    if (nrows == 1)
    {
        const size_t idx = arow_colind[0];
        const T val = arow_vals[0];

        auto browStart = bcolind.begin() + browptr[idx];
        const auto browEnd = bcolind.begin() + browptr[idx + 1];
        auto browStartVal = bvals.begin() + browptr[idx];

        while (browStart != browEnd)
        {
            *crow_colind++ = *browStart++;
            *crow_vals++ = val * (*browStartVal++);
        }

        return;
    }

    // Two rows, merge them
    if (nrows == 2)
    {
        const size_t row_colind1 = arow_colind[0];
        const size_t row_colind2 = arow_colind[1];
        const T row_val1 = arow_vals[0];
        const T row_val2 = arow_vals[1];

        std::span<const size_t> brow_colind1(bcolind.begin() + browptr[row_colind1], browptr[row_colind1 + 1] - browptr[row_colind1]);
        std::span<const size_t> brow_colind2(bcolind.begin() + browptr[row_colind2], browptr[row_colind2 + 1] - browptr[row_colind2]);
        std::span<const T> brow_vals1(bvals.begin() + browptr[row_colind1], browptr[row_colind1 + 1] - browptr[row_colind1]);
        std::span<const T> brow_vals2(bvals.begin() + browptr[row_colind2], browptr[row_colind2 + 1] - browptr[row_colind2]);
        mergeRows(row_val1, brow_colind1, brow_vals1, row_val2, brow_colind2, brow_vals2, crow_colind, crow_vals);

        return;
    }

    /**
     * Generic case (more than two rows).
     *
     * Merge rows by pairs, then merge the results together. When merging two rows, the result is always wider (or equal).
     * Merging by pairs allows to work with short rows as often as possible.
     */
    // Merge first two rows
    auto arow_colind_it = arow_colind.begin();
    size_t acol1 = *arow_colind_it++;
    size_t acol2 = *arow_colind_it++;

    auto arow_vals_it = arow_vals.begin();
    T aval1 = *arow_vals_it++;
    T aval2 = *arow_vals_it++;

    size_t *temp_col1 = crow_colind;
    T *temp_val1 = crow_vals;

    std::span<const size_t> brow_colind1(bcolind.begin() + browptr[acol1], browptr[acol1 + 1] - browptr[acol1]);
    std::span<const size_t> brow_colind2(bcolind.begin() + browptr[acol2], browptr[acol2 + 1] - browptr[acol2]);
    std::span<const T> brow_vals1(bvals.begin() + browptr[acol1], browptr[acol1 + 1] - browptr[acol1]);
    std::span<const T> brow_vals2(bvals.begin() + browptr[acol2], browptr[acol2 + 1] - browptr[acol2]);

    size_t c_numcol1 = mergeRows(aval1, brow_colind1, brow_vals1, aval2, brow_colind2, brow_vals2, temp_col1, temp_val1) - temp_col1;

    // Go by pairs
    while (arow_colind_it + 1 < arow_colind.end())
    {
        acol1 = *arow_colind_it++;
        acol2 = *arow_colind_it++;

        aval1 = *arow_vals_it++;
        aval2 = *arow_vals_it++;

        std::span<const size_t> brow_colindfirst(bcolind.begin() + browptr[acol1], browptr[acol1 + 1] - browptr[acol1]);
        std::span<const size_t> brow_colindsecond(bcolind.begin() + browptr[acol2], browptr[acol2 + 1] - browptr[acol2]);
        std::span<const T> brow_valsfirst(bvals.begin() + browptr[acol1], browptr[acol1 + 1] - browptr[acol1]);
        std::span<const T> brow_valssecond(bvals.begin() + browptr[acol2], browptr[acol2 + 1] - browptr[acol2]);

        size_t c_numcol2 = mergeRows(aval1, brow_colindfirst, brow_valsfirst, aval2, brow_colindsecond, brow_valssecond, temp_col2, temp_val2) - temp_col2;

        std::span<const size_t> temp_col1_span(temp_col1, c_numcol1);
        std::span<const size_t> temp_col2_span(temp_col2, c_numcol2);
        std::span<const T> temp_val1_span(temp_val1, c_numcol1);
        std::span<const T> temp_val2_span(temp_val2, c_numcol2);

        c_numcol1 = mergeRows(T(1), temp_col1_span, temp_val1_span, T(1), temp_col2_span, temp_val2_span, temp_col3, temp_val3) - temp_col3;

        std::swap(temp_col3, temp_col1);
        std::swap(temp_val3, temp_val1);
    }

    // Merge the tail
    if (arow_colind_it < arow_colind.end())
    {
        acol2 = *arow_colind_it++;
        aval2 = *arow_vals_it++;
        std::span<const size_t> brow_tail(bcolind.begin() + browptr[acol2], browptr[acol2 + 1] - browptr[acol2]);
        std::span<const T> brow_vals_tail(bvals.begin() + browptr[acol2], browptr[acol2 + 1] - browptr[acol2]);

        std::span<const size_t> temp_col1_span(temp_col1, c_numcol1);
        std::span<const T> temp_val1_span(temp_val1, c_numcol1);
        c_numcol1 = mergeRows(T(1), temp_col1_span, temp_val1_span, aval2, brow_tail, brow_vals_tail, temp_col3, temp_val3) - temp_col3;

        std::swap(temp_col3, temp_col1);
        std::swap(temp_val3, temp_val1);
    }

    if (temp_col1 != crow_colind)
    {
        std::copy(temp_col1, temp_col1 + c_numcol1, crow_colind);
        std::copy(temp_val1, temp_val1 + c_numcol1, crow_vals);
    }
}

template <typename T>
void spgemm_rmerge(const CSRMatrix<T> &A, const CSRMatrix<T> &B, CSRMatrix<T> &C)
{
    std::atomic_size_t maxRowWidth = 0;

    auto countRowWidth = [&](size_t start, size_t end) {
        size_t threadMax = 0;
        for (size_t i = start; i < end; ++i) {
            const size_t rowStart = (*(A.row_pointers))[i];
            const size_t rowEnd = (*(A.row_pointers))[i + 1];
            size_t rowWidth = 0;

            for (size_t j = rowStart; j < rowEnd; ++j) {
                size_t colIdx = (*(A.col_indices))[j];
                rowWidth += (*(B.row_pointers))[colIdx + 1] - (*(B.row_pointers))[colIdx];
            }

            threadMax = std::max(threadMax, rowWidth);
        }
        
        if (threadMax > maxRowWidth) {
            maxRowWidth = threadMax;
        }
    };

    launchThreads(A.row_num, countRowWidth);

    // Temporary row of C for each thread. For now the number of threads is 1
    std::vector<std::vector<size_t>> tempCol(num_threads);
    std::vector<std::vector<T>> tempVal(num_threads);
    for (int i = 0; i < num_threads; ++i) {
        tempCol[i].resize(3 * maxRowWidth);
        tempVal[i].resize(2 * maxRowWidth);
    }

    C.row_num = A.row_num;
    C.col_num = B.col_num;
    C.row_pointers = new std::vector<size_t>(C.row_num + 1, 0);

    auto fillRowPointers = [&](size_t start, size_t end, int thread_id) {
        for (size_t i = start; i < end; ++i) {
            const size_t rowStart = (*(A.row_pointers))[i];
            const size_t rowEnd = (*(A.row_pointers))[i + 1];
            std::span<const size_t> arow(A.col_indices->begin() + rowStart, rowEnd - rowStart);
            (*(C.row_pointers))[i + 1] = prodRowWidth(arow, *(B.row_pointers), *(B.col_indices), tempCol[thread_id].data(), tempCol[thread_id].data() + maxRowWidth, tempCol[thread_id].data() + 2 * maxRowWidth);
        }
    };

    launchThreadsWithID(C.row_num, fillRowPointers);

    C.nnz = C.scanRowSize();
    C.col_indices = new std::vector<size_t>(C.nnz, 0);
    C.vals = new std::vector<T>(C.nnz, T());

    auto computeColIndices = [&](size_t start, size_t end, int thread_id) {
        for (size_t i = start; i < end; ++i) {
            const size_t rowStart = (*(A.row_pointers))[i];
            const size_t rowEnd = (*(A.row_pointers))[i + 1];
            std::span<const size_t> arow_colind(A.col_indices->begin() + rowStart, rowEnd - rowStart);
            std::span<const T> arow_vals(A.vals->begin() + rowStart, rowEnd - rowStart);

            size_t *crow_colind = C.col_indices->data() + (*(C.row_pointers))[i];
            T *crow_vals = C.vals->data() + (*(C.row_pointers))[i];

            prodRow(arow_colind, arow_vals, *(B.row_pointers), *(B.col_indices), *(B.vals),
                    crow_colind, crow_vals,
                    tempCol[thread_id].data(), tempVal[thread_id].data(),
                    tempCol[thread_id].data() + maxRowWidth, tempVal[thread_id].data() + maxRowWidth);
        }
    };

    launchThreadsWithID(C.row_num, computeColIndices);
}