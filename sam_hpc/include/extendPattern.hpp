#pragma once

#include "CSRMatrix.hpp"

void saad_pattern_extension(CSRMatrix<int> &prevPattern, CSRMatrix<int> &nextPattern)
{
    nextPattern.row_num = prevPattern.col_num;
    nextPattern.col_num = prevPattern.row_num;
    nextPattern.row_pointers = new std::vector<size_t>(prevPattern.row_num + 1);

    // Count the number of nonzero elements in each row of nextPattern
    std::vector<int> marker(nextPattern.col_num, -1);
    for (size_t ia = 0; ia < prevPattern.row_num; ++ia)
    {
        int count = 0;
        const size_t rowBegA = (*(prevPattern.row_pointers))[ia];
        const size_t rowEndA = (*(prevPattern.row_pointers))[ia + 1];
        for (size_t ja = rowBegA; ja < rowEndA; ++ja)
        {
            const size_t colIdxA = (*(prevPattern.col_indices))[ja];
            const size_t rowBegB = (*(prevPattern.row_pointers))[colIdxA];
            const size_t rowEndB = (*(prevPattern.row_pointers))[colIdxA + 1];

            for (size_t jb = rowBegB; jb < rowEndB; ++jb)
            {
                const size_t colIdxB = (*(prevPattern.col_indices))[jb];
                if (marker[colIdxB] != static_cast<int>(ia))
                {
                    marker[colIdxB] = static_cast<int>(ia);
                    ++count;
                }
            }
        }
        (*(nextPattern.row_pointers))[ia + 1] = count;
    }

    nextPattern.nnz = nextPattern.scanRowSize();
    nextPattern.col_indices = new std::vector<size_t>(nextPattern.nnz, 0);
    nextPattern.vals = new std::vector<int>(nextPattern.nnz, 1);

    // Compute the column indices
    marker.assign(nextPattern.col_num, -1);
    for (size_t ia = 0; ia < nextPattern.row_num; ++ia)
    {
        const size_t rowBeg = (*(nextPattern.row_pointers))[ia];
        size_t rowEnd = rowBeg;
        for (size_t ja = (*(prevPattern.row_pointers))[ia], ea = (*(prevPattern.row_pointers))[ia + 1]; ja < ea; ++ja)
        {
            size_t colIdxA = (*(prevPattern.col_indices))[ja];

            for (size_t jb = (*(prevPattern.row_pointers))[colIdxA], eb = (*(prevPattern.row_pointers))[colIdxA + 1]; jb < eb; ++jb)
            {
                size_t colIdxB = (*(prevPattern.col_indices))[jb];
                if (marker[colIdxB] < static_cast<int>(rowBeg))
                {
                    marker[colIdxB] = static_cast<int>(rowEnd);
                    (*(nextPattern.col_indices))[rowEnd] = colIdxB;
                    ++rowEnd;
                }
            }
        }

        std::sort(nextPattern.col_indices->begin() + rowBeg, nextPattern.col_indices->begin() + rowEnd);
    }
}

template <bool needOut>
size_t* mergeRows(std::span<const size_t> row1, std::span<const size_t> row2, size_t* result)
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
    return size_t(); // Control should never reach here
}

size_t* mergeRowsComputeIndices(std::span<const size_t> row_colind1, std::span<const size_t> row_colind2, size_t* crow_colind)
{
    auto row_colind1_it = row_colind1.begin();
    auto row_colind2_it = row_colind2.begin();

    while (row_colind1_it != row_colind1.end() && row_colind2_it != row_colind2.end())
    {
        const size_t col1 = *row_colind1_it;
        const size_t col2 = *row_colind2_it;

        if (col1 < col2)
        {
            ++row_colind1_it;
            *crow_colind = col1;
        }
        else if (col1 == col2)
        {
            ++row_colind1_it;
            ++row_colind2_it;

            *crow_colind = col1;
        }
        else
        {
            ++row_colind2_it;
            *crow_colind = col2;
        }

        ++crow_colind;
    }

    while (row_colind1_it < row_colind1.end())
    {
        *crow_colind++ = *row_colind1_it++;
    }

    while (row_colind2_it < row_colind2.end())
    {
        *crow_colind++ = *row_colind2_it++;
    }

    return crow_colind;
}

size_t prodRowWidth(std::span<const size_t> arow, const std::vector<size_t> &browptr, const std::vector<size_t> &bcolind,
                    size_t* temp_col1, size_t* temp_col2, size_t* temp_col3)
{
    const size_t nrows = arow.size();

    // No rows merge, nothing to do
    if (nrows == 0)
        return 0;

    // Single row, just copy it to output
    if (nrows == 1)
        return browptr[arow[0] + 1] - browptr[arow[0]];

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

        if (arow_it == arow.end())
        {
            std::span<const size_t> temp_col1_subspan(temp_col1, ncols1);
            std::span<const size_t> temp_col2_subspan(temp_col2, ncols2);
            return mergeRows<false>(temp_col1_subspan, temp_col2_subspan, temp_col3) - temp_col3;
        }
        else
        {
            std::span<const size_t> temp_col1_subspan(temp_col1, ncols1);
            std::span<const size_t> temp_col2_subspan(temp_col2, ncols2);
            ncols1 = mergeRows<true>(temp_col1_subspan, temp_col2_subspan, temp_col3) - temp_col3;
            std::swap(temp_col1, temp_col3);
        }
    }

    // Merge the tail
    const auto tail = *arow_it++;
    std::span<const size_t> browtail(bcolind.begin() + browptr[tail], browptr[tail + 1] - browptr[tail]);
    std::span<const size_t> temp_col1_subspan(temp_col1, ncols1);
    return mergeRows<false>(temp_col1_subspan, browtail, temp_col2) - temp_col2;
}

void prodRow(std::span<const size_t> arow_colind, const std::vector<size_t> &browptr, const std::vector<size_t> &bcolind,
             size_t* crow_colind, size_t* temp_col2, size_t* temp_col3)
{
    const size_t nrows = arow_colind.size();

    // No rows to merge, nothing to do
    if (nrows == 0)
        return;

    // Single row, just copy it to output
    if (nrows == 1)
    {
        const size_t idx = arow_colind[0];

        auto browStart = bcolind.begin() + browptr[idx];
        const auto browEnd = bcolind.begin() + browptr[idx + 1];

        while (browStart != browEnd)
        {
            *crow_colind++ = *browStart++;
        }

        return;
    }

    // Two rows, merge them
    if (nrows == 2)
    {
        const size_t row_colind1 = arow_colind[0];
        const size_t row_colind2 = arow_colind[1];

        std::span<const size_t> brow_colind1(bcolind.begin() + browptr[row_colind1], browptr[row_colind1 + 1] - browptr[row_colind1]);
        std::span<const size_t> brow_colind2(bcolind.begin() + browptr[row_colind2], browptr[row_colind2 + 1] - browptr[row_colind2]);
        mergeRowsComputeIndices(brow_colind1, brow_colind2, crow_colind);

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

    size_t *temp_col1 = crow_colind;

    std::span<const size_t> brow_colind1(bcolind.begin() + browptr[acol1], browptr[acol1 + 1] - browptr[acol1]);
    std::span<const size_t> brow_colind2(bcolind.begin() + browptr[acol2], browptr[acol2 + 1] - browptr[acol2]);

    size_t c_numcol1 = mergeRowsComputeIndices(brow_colind1, brow_colind2, temp_col1) - temp_col1;

    // Go by pairs
    while (arow_colind_it + 1 < arow_colind.end())
    {
        acol1 = *arow_colind_it++;
        acol2 = *arow_colind_it++;

        std::span<const size_t> brow_colindfirst(bcolind.begin() + browptr[acol1], browptr[acol1 + 1] - browptr[acol1]);
        std::span<const size_t> brow_colindsecond(bcolind.begin() + browptr[acol2], browptr[acol2 + 1] - browptr[acol2]);

        size_t c_numcol2 = mergeRowsComputeIndices(brow_colindfirst, brow_colindsecond, temp_col2) - temp_col2;

        std::span<const size_t> temp_col1_subspan(temp_col1, c_numcol1);
        std::span<const size_t> temp_col2_subspan(temp_col2, c_numcol2);
        c_numcol1 = mergeRowsComputeIndices(temp_col1_subspan, temp_col2_subspan, temp_col3) - temp_col3;

        std::swap(temp_col3, temp_col1);
    }

    // Merge the tail
    if (arow_colind_it < arow_colind.end())
    {
        acol2 = *arow_colind_it++;
        std::span<const size_t> brow_tail(bcolind.begin() + browptr[acol2], browptr[acol2 + 1] - browptr[acol2]);

        std::span<size_t> temp_col1_subspan(temp_col1, c_numcol1);
        c_numcol1 = mergeRowsComputeIndices(temp_col1_subspan, brow_tail, temp_col3) - temp_col3;

        std::swap(temp_col3, temp_col1);
    }

    if (temp_col1 != crow_colind)
    {
        std::copy(temp_col1, temp_col1 + c_numcol1, crow_colind);
    }
}

void rmerge_pattern_extension(CSRMatrix<int> &prevPattern, CSRMatrix<int> &nextPattern)
{
    size_t maxRowWidth = 0;

    //@todo: write the parallel version
    for (size_t i = 0; i < prevPattern.row_num; ++i)
    {
        const size_t rowStart = (*(prevPattern.row_pointers))[i];
        const size_t rowEnd = (*(prevPattern.row_pointers))[i + 1];
        size_t rowWidth = 0;

        for (size_t j = rowStart; j < rowEnd; ++j)
        {
            size_t colIdx = (*(prevPattern.col_indices))[j];
            rowWidth += (*(prevPattern.row_pointers))[colIdx + 1] - (*(prevPattern.row_pointers))[colIdx];
        }

        maxRowWidth = std::max(maxRowWidth, rowWidth);
    }

    // Temporary row of C for each thread. For now the number of threads is 1
    std::vector<std::vector<size_t>> tempCol(1, std::vector<size_t>(3 * maxRowWidth, 0));

    nextPattern.row_num = prevPattern.row_num;
    nextPattern.col_num = prevPattern.col_num;
    nextPattern.row_pointers = new std::vector<size_t>(nextPattern.row_num + 1, 0);

    for (size_t i = 0; i < nextPattern.row_num; ++i)
    {
        const size_t rowStart = (*(prevPattern.row_pointers))[i];
        const size_t rowEnd = (*(prevPattern.row_pointers))[i + 1];
        std::span<const size_t> arow(prevPattern.col_indices->data() + rowStart, rowEnd - rowStart);
        (*(nextPattern.row_pointers))[i + 1] = prodRowWidth(arow, *(prevPattern.row_pointers), *(prevPattern.col_indices),
                                                            tempCol[0].data(), tempCol[0].data() + maxRowWidth, tempCol[0].data() + 2 * maxRowWidth);
    }

    nextPattern.nnz = nextPattern.scanRowSize();
    nextPattern.col_indices = new std::vector<size_t>(nextPattern.nnz, 0);
    nextPattern.vals = new std::vector<int>(nextPattern.nnz, 1);

    for (size_t i = 0; i < nextPattern.row_num; ++i)
    {
        const size_t rowStart = (*(prevPattern.row_pointers))[i];
        const size_t rowEnd = (*(prevPattern.row_pointers))[i + 1];
        std::span<const size_t> arow_colind(prevPattern.col_indices->begin() + rowStart, rowEnd - rowStart);

        size_t *crow_colind = nextPattern.col_indices->data() + (*(nextPattern.row_pointers))[i];

        prodRow(arow_colind, *(prevPattern.row_pointers), *(prevPattern.col_indices), crow_colind, tempCol[0].data(), tempCol[0].data() + maxRowWidth);
    }
}

void extend_pattern(CSRMatrix<int> &S, const int level)
{
    std::vector<CSRMatrix<int>> finalPatterns(level - 1);

    for (int l = 0; l < level - 1; ++l)
    {
        if (l == 0)
        {
            // saad_pattern_extension(S, finalPatterns[l]);
            rmerge_pattern_extension(S, finalPatterns[l]);
        }
        else if (l >= 1)
        {
            // saad_pattern_extension(finalPatterns[l - 1], finalPatterns[l]);
            rmerge_pattern_extension(finalPatterns[l - 1], finalPatterns[l]);
        }
    }

    swap(S, finalPatterns.back());
}