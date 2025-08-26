#pragma once

#include "ilu_solve.hpp"

/// ILU(0) smoother.
template <typename T>
struct ilu0 {
    /// Relaxation parameters.
    struct params {
        /// Damping factor.
        T damping;

        /// Parameters for sparse triangular system solver
        // typename ilu_solve_iterative<T>::params solve;
        typename ilu_solve_direct<T>::params solve;
        
        params() : damping(1) {}
    } prm;

    ilu0(const CSRMatrix<T> &A, const params &prm) : prm(prm) {
        const size_t n = A.m_rows;

        size_t Lnz = 0, Unz = 0;

        for(size_t i = 0; i < n; ++i) {
            size_t row_beg = A.m_row_pointers[i];
            size_t row_end = A.m_row_pointers[i + 1];

            for(size_t j = row_beg; j < row_end; ++j) {
                size_t c = A.m_col_indices[j];
                if (c < i)
                    ++Lnz;
                else if (c > i)
                    ++Unz;
            }
        }

        std::shared_ptr<CSRMatrix<T>> L = std::make_shared<CSRMatrix<T>>();
        std::shared_ptr<CSRMatrix<T>> U = std::make_shared<CSRMatrix<T>>();

        L->m_rows = n;
        L->m_cols = n;
        L->m_nnz = Lnz;
        L->m_row_pointers.resize(L->m_rows + 1);
        L->m_col_indices.resize(L->m_nnz);
        L->m_vals.resize(L->m_nnz, T(0));

        U->m_rows = n;
        U->m_cols = n;
        U->m_nnz = Unz;
        U->m_row_pointers.resize(U->m_rows + 1);
        U->m_col_indices.resize(U->m_nnz);
        U->m_vals.resize(U->m_nnz, T(0));

        size_t Lhead = 0;
        size_t Uhead = 0;

        std::shared_ptr<std::vector<T>> D = std::make_shared<std::vector<T>>(n, T(0));

        std::vector<T*> work(n, NULL);

        for(size_t i = 0; i < n; ++i) {
            size_t row_beg = A.m_row_pointers[i];
            size_t row_end = A.m_row_pointers[i + 1];

            for(size_t j = row_beg; j < row_end; ++j) {
                size_t  c = A.m_col_indices[j];
                T v = A.m_vals[j];

                if (c < i) {
                    L->m_col_indices[Lhead] = c;
                    L->m_vals[Lhead] = v;
                    work[c] = L->m_vals.data() + Lhead;
                    ++Lhead;
                } else if (c == i) {
                    (*D)[i] = v + 1e-10;
                    work[c] = &(*D)[i];
                } else {
                    U->m_col_indices[Uhead] = c;
                    U->m_vals[Uhead] = v;
                    work[c] = U->m_vals.data() + Uhead;
                    ++Uhead;
                }
            }

            L->m_row_pointers[i + 1] = Lhead;
            U->m_row_pointers[i + 1] = Uhead;

            for(size_t j = row_beg; j < row_end; ++j) {
                size_t c = A.m_col_indices[j];

                // Exit if diagonal is reached
                if (c >= i) {
                    assert(c == i && "No diagonal value in system matrix");
                    // assert((*D)[i] == T(0) && "Zero pivot in ILU");
                    if ((*D)[i] == T(0))
                        (*D)[i] = 1e-10;

                    (*D)[i] = 1 / (*D)[i];
                    break;
                }

                // Compute the multiplier for jrow
                T tl = (*work[c]) * (*D)[c];
                *work[c] = tl;

                // Perform linear combination
                for(size_t k = U->m_row_pointers[c]; k < U->m_row_pointers[c + 1]; ++k) {
                    T *w = work[U->m_col_indices[k]];
                    if (w) *w -= tl * U->m_vals[k];
                }
            }

            // Get rid of zeros in the factors
            Lhead = L->m_row_pointers[i];
            Uhead = U->m_row_pointers[i];

            for(size_t j = Lhead, e = L->m_row_pointers[i + 1]; j < e; ++j) {
                T v = L->m_vals[j];
                if (v != T(0)) {
                    L->m_col_indices[Lhead] = L->m_col_indices[j];
                    L->m_vals[Lhead] = v;
                    ++Lhead;
                }
            }

            for(size_t j = Uhead, e = U->m_row_pointers[i + 1]; j < e; ++j) {
                T v = U->m_vals[j];
                if (v != T(0)) {
                    U->m_col_indices[Uhead] = U->m_col_indices[j];
                    U->m_vals[Uhead] = v;
                    ++Uhead;
                }
            }
            L->m_row_pointers[i + 1] = Lhead;
            U->m_row_pointers[i + 1] = Uhead;

            // Refresh work
            for(size_t j = row_beg; j < row_end; ++j)
                work[A.m_col_indices[j]] = NULL;
        }

        L->m_nnz = Lhead;
        U->m_nnz = Uhead;

        // ilu = std::make_shared<ilu_solve_iterative<T>>(L, U, D, prm.solve);
        ilu = std::make_shared<ilu_solve_direct<T>>(L, U, D, prm.solve);
    }

    void apply_pre(const CSRMatrix<T> &A, const std::vector<T> &rhs,
                   std::vector<T> &x, std::vector<T> &tmp) const
    {
        residual(rhs, A, x, tmp);
        ilu->solve(tmp);
        axpby(prm.damping, tmp, T(1), x);
    }

    void apply_post(const CSRMatrix<T> &A, const std::vector<T> &rhs,
                    std::vector<T> &x, std::vector<T> &tmp) const
    {
        residual(rhs, A, x, tmp);
        ilu->solve(tmp);
        axpby(prm.damping, tmp, T(1), x);
    }

    void apply(const CSRMatrix<T> &A, const std::vector<T> &rhs, std::vector<T> &x) const {
        x = rhs;
        ilu->solve(x);
    }

    // Todo: bytes method

    private:
        // std::shared_ptr<ilu_solve_iterative<T>> ilu;
        std::shared_ptr<ilu_solve_direct<T>> ilu;
};