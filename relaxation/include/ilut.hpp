#pragma once

#include "ilu_solve.hpp"

#include <queue>

/// ILUT(p, tau) smoother.
template <typename T>
struct ilut {
    /// Relaxation parameters.
    struct params {
        /// Fill factor.
        T p;

        /// Minimum magnitude of non-zero elements relative to the current row norm.
        T tau;

        /// Damping factor.
        T damping;

        /// Parameters for sparse triangular system solver
        // typename ilu_solve_iterative<T>::params solve;
        typename ilu_solve_direct<T>::params solve;

        params() : p(2), tau(1e-2f), damping(1) {}
    } prm;

    ilut(const CSRMatrix<T> &A, const params &prm) : prm(prm) {
        const size_t n = A.m_rows;

        size_t Lnz = 0, Unz = 0;

        for(size_t i = 0; i < n; ++i) {
            size_t row_beg = A.m_row_pointers[i];
            size_t row_end = A.m_row_pointers[i + 1];

            int lenL = 0, lenU = 0;
            for(size_t j = row_beg; j < row_end; ++j) {
                size_t c = A.m_col_indices[j];
                if (c < i)
                    ++lenL;
                else if (c > i)
                    ++lenU;
            }

            Lnz += static_cast<size_t>(lenL * prm.p);
            Unz += static_cast<size_t>(lenU * prm.p);
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

        std::shared_ptr<std::vector<T>> D = std::make_shared<std::vector<T>>(n, T(0));

        sparse_vector w(n);

        for(size_t i = 0, Lhead = 0, Uhead = 0; i < n; ++i) {
            w.dia = i;

            int lenL = 0;
            int lenU = 0;

            T tol = T(0);

            for (auto a = A.rowBegin(i); a; ++a) {
                w[a.col()] = a.value();
                tol += std::abs(a.value());

                if (a.col() < i) ++lenL;
                if (a.col() > i) ++lenU;
            }
            tol *= prm.tau / (lenL + lenU);

            while(!w.q.empty()) {
                ptrdiff_t k = w.next_nonzero();
                w[k] = w[k] * (*D)[k];
                T wk = w[k];

                if (std::abs(wk) > tol) {
                    for(size_t j = U->m_row_pointers[k]; j < U->m_row_pointers[k + 1]; ++j)
                        w[U->m_col_indices[j]] -= wk * U->m_vals[j];
                }
            }

            w.move_to(static_cast<int>(lenL * prm.p), static_cast<int>(lenU * prm.p), tol, Lhead, *L, Uhead, *U, *D);

            L->m_row_pointers[i + 1] = Lhead;
            U->m_row_pointers[i + 1] = Uhead;
        }

        L->m_nnz = L->m_row_pointers[n];
        U->m_nnz = U->m_row_pointers[n];

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

    // TODO: bytes method

    private:
        // std::shared_ptr<ilu_solve_iterative<T>> ilu;
        std::shared_ptr<ilu_solve_direct<T>> ilu;

        struct sparse_vector {
            struct nonzero {
                ptrdiff_t  col;
                T val;

                nonzero() : col(-1) {}

                nonzero(ptrdiff_t col, const T &val = T(0))
                    : col(col), val(val) {}
            };

            struct comp_indices {
                const std::vector<nonzero> &nz;

                comp_indices(const std::vector<nonzero> &nz) : nz(nz) {}

                bool operator()(int a, int b) const {
                    return nz[a].col > nz[b].col;
                }
            };

            typedef std::priority_queue<int, std::vector<int>, comp_indices> priority_queue;

            std::vector<nonzero>   nz;
            std::vector<ptrdiff_t> idx;
            priority_queue q;

            ptrdiff_t dia;

            sparse_vector(size_t n) : idx(n, -1), q(comp_indices(nz)), dia(0) {
                nz.reserve(16);
            }

            T operator[](ptrdiff_t i) const {
                if (idx[i] >= 0) return nz[idx[i]].val;
                return T(0);
            }

            T& operator[](ptrdiff_t i) {
                if (idx[i] == -1) {
                    int p = nz.size();
                    idx[i] = p;
                    nz.push_back(nonzero(i));
                    if (i < dia) q.push(p);
                }
                return nz[idx[i]].val;
            }

            typename std::vector<nonzero>::iterator begin() {
                return nz.begin();
            }

            typename std::vector<nonzero>::iterator end() {
                return nz.end();
            }

            ptrdiff_t next_nonzero() {
                int p = q.top();
                q.pop();
                return nz[p].col;
            }

            struct higher_than {
                T tol;
                ptrdiff_t dia;

                higher_than(T tol, ptrdiff_t dia)
                    : tol(tol), dia(dia) {}

                bool operator()(const nonzero &v) const {
                    return v.col == dia || std::abs(v.val) > tol;
                }
            };

            struct L_first {
                ptrdiff_t dia;

                L_first(ptrdiff_t dia) : dia(dia) {}

                bool operator()(const nonzero &v) const {
                    return v.col < dia;
                }
            };

            struct by_abs_val {
                ptrdiff_t dia;

                by_abs_val(ptrdiff_t dia) : dia(dia) {}

                bool operator()(const nonzero &a, const nonzero &b) const {
                    if (a.col == dia) return true;
                    if (b.col == dia) return false;

                    return std::abs(a.val) > std::abs(b.val);
                }
            };

            struct by_col {
                bool operator()(const nonzero &a, const nonzero &b) const {
                    return a.col < b.col;
                }
            };

            void move_to(int lp, int up, T tol, size_t &Lhead, CSRMatrix<T> &L, size_t &Uhead, CSRMatrix<T> &U, std::vector<T> &D) {
                typedef typename std::vector<nonzero>::iterator ptr;

                ptr b = nz.begin();
                ptr e = nz.end();

                // Move zeros to back:
                e = std::partition(b, e, higher_than(tol, dia));

                // Split L and U:
                ptr m = std::partition(b, e, L_first(dia));

                // Get largest p elements in L and U.
                ptr lend = std::min(b + lp, m);
                ptr uend = std::min(m + up, e);

                if (lend != m) std::nth_element(b, lend, m, by_abs_val(dia));
                if (uend != e) std::nth_element(m, uend, e, by_abs_val(dia));

                // Sort entries by column number
                std::sort(b, lend, by_col());
                std::sort(m, uend, by_col());

                // copy L to the output matrix.
                for(ptr a = b; a != lend; ++a) {
                    L.m_col_indices[Lhead] = a->col;
                    L.m_vals[Lhead] = a->val;

                    ++Lhead;
                }

                // Store inverted diagonal.
                D[dia] = 1 / m->val;

                if (m != uend) {
                    ++m;

                    // copy U to the output matrix.
                    for(ptr a = m; a != uend; ++a) {
                        U.m_col_indices[Uhead] = a->col;
                        U.m_vals[Uhead] = a->val;

                        ++Uhead;
                    }
                }

                for(const nonzero &e : nz) idx[e.col] = -1;
                nz.clear();
            }
        };
};