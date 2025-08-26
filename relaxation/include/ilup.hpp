#pragma once

#include "ilu_solve.hpp"
#include "ilu0.hpp"

template <typename T>
std::shared_ptr<CSRMatrix<T>> symb_product(const CSRMatrix<T> &A, const CSRMatrix<T> &B) {
    std::shared_ptr<CSRMatrix<T>> C = std::make_shared<CSRMatrix<T>>();

    C->m_rows = A.m_rows;
    C->m_cols = A.m_cols;
    C->m_row_pointers.resize(C->m_rows);

    tbb::enumerable_thread_specific<std::vector<ptrdiff_t>> local_markers_1(B.m_cols, -1);
    tbb::parallel_for(tbb::blocked_range<size_t>(0, A.m_rows), [&](const tbb::blocked_range<size_t>& r) {
        std::vector<ptrdiff_t>& marker = local_markers_1.local();
        for (size_t ia = r.begin(); ia < r.end(); ++ia) {
            int count = 0;
            for (size_t ja = A.m_row_pointers[ia], ea = A.m_row_pointers[ia + 1]; ja < ea; ++ja) {
                const size_t colIdxA = A.m_col_indices[ja];

                for (size_t jb = B.m_row_pointers[colIdxA], eb = B.m_row_pointers[colIdxA + 1]; jb < eb; ++jb) {
                    const size_t colIdxB = B.m_col_indices[jb];
                    if (marker[colIdxB] != static_cast<ptrdiff_t>(ia)) {
                        marker[colIdxB] = static_cast<ptrdiff_t>(ia);
                        ++count;
                    }
                }
            }
            C->m_row_pointers[ia + 1] = count;
        }
    });

    C->m_nnz = C->scanRowSize();
    C->m_col_indices.resize(C->m_nnz, 0);
    C->m_vals.resize(C->m_nnz, T(0));

    // Compute the column indices
    local_markers_1.clear();
    // tbb::enumerable_thread_specific<std::vector<int>> local_markers_2(nextPattern.m_cols, -1);
    tbb::parallel_for(tbb::blocked_range<size_t>(0, C->m_rows), [&](const tbb::blocked_range<size_t>& r) {
        std::vector<ptrdiff_t>& marker = local_markers_1.local();
        for (size_t ia = r.begin(); ia < r.end(); ++ia) {
            const size_t rowBeg = C->m_row_pointers[ia];
            size_t rowEnd = rowBeg;

            for (size_t ja = A.m_row_pointers[ia], ea = A.m_row_pointers[ia + 1]; ja < ea; ++ja) {
                size_t colIdxA = A.m_col_indices[ja];

                for (size_t jb = B.m_row_pointers[colIdxA], eb = B.m_row_pointers[colIdxA + 1]; jb < eb; ++jb) {
                    size_t colIdxB = B.m_col_indices[jb];
                    if (marker[colIdxB] < static_cast<ptrdiff_t>(rowBeg)) {
                        marker[colIdxB] = static_cast<ptrdiff_t>(rowEnd);
                        C->m_col_indices[rowEnd] = colIdxB;
                        ++rowEnd;
                    }
                }
            }

            for (size_t k = rowBeg; k < rowEnd; ++k) {
                marker[C->m_col_indices[k]] = -1;
            }

            std::sort(C->m_col_indices.begin() + rowBeg, C->m_col_indices.begin() + rowEnd);
        }
    });

    return C;
}

/// ILU(k) smoother.
template <class T>
struct ilup {
    typedef ilu0<T> Base;

    /// Relaxation parameters.
    struct params : Base::params {
        typedef typename Base::params BasePrm;

        /// Level of fill-in.
        int k;

        params() : k(1) {}
    } prm;

    ilup(const CSRMatrix<T> &A, const params &prm) : prm(prm) {
        if (prm.k == 0) {
            base = std::make_shared<Base>(A, prm);
        } else {
            std::shared_ptr<CSRMatrix<T>> P = symb_product(A, A);
            for(int k = 1; k < prm.k; ++k) {
                P = symb_product(*P, A);
            }

            size_t n = A.m_rows;

            tbb::parallel_for(tbb::blocked_range<size_t>(0, n), [&](const tbb::blocked_range<size_t> &r) {
                for(size_t i = r.begin(); i < r.end(); ++i) {
                    size_t p_beg = P->m_row_pointers[i];
                    size_t p_end = P->m_row_pointers[i + 1];
                    size_t a_beg = A.m_row_pointers[i];
                    size_t a_end = A.m_row_pointers[i + 1];

                    // std::fill(P->m_vals.data() + p_beg, P->m_vals.data() + p_end, T(0));

                    for(ptrdiff_t ja = a_beg, ea = a_end, jp = p_beg, ep = p_end; ja < ea; ++ja) {
                        ptrdiff_t ca = A.m_col_indices[ja];
                        while(jp < ep && P->m_col_indices[jp] < ca) ++jp;
                        if (P->m_col_indices[jp] == ca) P->m_vals[jp] = A.m_vals[ja];
                    }
                }
            });

            base = std::make_shared<Base>(*P, prm);
        }
    }

    void apply_pre(const CSRMatrix<T> &A, const std::vector<T> &rhs,
                   std::vector<T> &x, std::vector<T> &tmp) const
    {
        base->apply_pre(A, rhs, x, tmp);
    }

    void apply_post(const CSRMatrix<T> &A, const std::vector<T> &rhs,
                   std::vector<T> &x, std::vector<T> &tmp) const
    {
        base->apply_post(A, rhs, x, tmp);
    }

    void apply(const CSRMatrix<T> &A, const std::vector<T> &rhs, std::vector<T> &x) const {
        base->apply(A, rhs, x);
    }

    // TODO: bytes method

    private:
        std::shared_ptr<Base> base;
};