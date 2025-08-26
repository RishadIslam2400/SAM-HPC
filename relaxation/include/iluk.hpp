#pragma once

#include "ilu_solve.hpp"

#include <queue>


/// ILU(k) smoother.
template <typename T>
struct iluk {
    /// Relaxation parameters.
    struct params {
        /// Level of fill-in.
        int k;

        /// Damping factor.
        T damping;

        /// Parameters for sparse triangular system solver
        // typename ilu_solve_iterative<T>::params solve;
        typename ilu_solve_direct<T>::params solve;

        params() : k(1), damping(1) {}
    } prm;

    iluk(const CSRMatrix<T> &A, const params &prm) : prm(prm) {
        const size_t n = A.m_rows;

        size_t Anz = A.m_nnz;

        std::vector<size_t> Lptr; 
        Lptr.reserve(n + 1); 
        Lptr.push_back(0);
        std::vector<size_t> Lcol; 
        Lcol.reserve(Anz / 3);
        std::vector<T> Lval; 
        Lval.reserve(Anz / 3);

        std::vector<size_t> Uptr; 
        Uptr.reserve(n + 1);
        Uptr.push_back(0);
        std::vector<size_t> Ucol;
        Ucol.reserve(Anz / 3);
        std::vector<T> Uval; 
        Uval.reserve(Anz / 3);

        std::vector<int> Ulev; 
        Ulev.reserve(Anz / 3);

        std::shared_ptr<std::vector<T>> D = std::make_shared<std::vector<T>>(n);

        sparse_vector w(n, prm.k);

        for(size_t i = 0; i < n; ++i) {
            w.reset(i);

            for(auto a = A.rowBegin(i); a; ++a) {
                w.add(a.col(), a.value(), 0);
            }

            while(!w.q.empty()) {
                nonzero &a = w.next_nonzero();
                a.val = a.val * (*D)[a.col];

                for(size_t j = Uptr[a.col], e = Uptr[a.col+1]; j < e; ++j) {
                    int lev = std::max(a.lev, Ulev[j]) + 1;
                    w.add(Ucol[j], -a.val * Uval[j], lev);
                }
            }

            w.sort();

            for(const nonzero &e : w.nz) {
                if (e.col < i) {
                    Lcol.push_back(e.col);
                    Lval.push_back(e.val);
                } else if (e.col == i) {
                    (*D)[i] = 1 / e.val;
                } else {
                    Ucol.push_back(e.col);
                    Uval.push_back(e.val);
                    Ulev.push_back(e.lev);
                }
            }

            Lptr.push_back(Lcol.size());
            Uptr.push_back(Ucol.size());
        }

        /* ilu = std::make_shared<ilu_solve_iterative<T>>(
                std::make_shared<CSRMatrix<T>>(n, n, Lptr, Lcol, Lval),
                std::make_shared<CSRMatrix<T>>(n, n, Uptr, Ucol, Uval),
                D, prm.solve); */

        ilu = std::make_shared<ilu_solve_direct<T>>(
                std::make_shared<CSRMatrix<T>>(n, n, Lval, Lptr, Lcol),
                std::make_shared<CSRMatrix<T>>(n, n, Uval, Uptr, Ucol),
                D, prm.solve);
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

    void apply(const CSRMatrix<T>&, const std::vector<T> &rhs, std::vector<T> &x) const {
        x = rhs;
        ilu->solve(x);
    }

    // TODO: bytes method

private:
    // std::shared_ptr<ilu_solve_iterative<T>> ilu;
    std::shared_ptr<ilu_solve_direct<T>> ilu;

    struct nonzero {
        ptrdiff_t  col;
        T val;
        int lev;

        nonzero() : col(-1) {}

        nonzero(ptrdiff_t col, const T &val, int lev)
            : col(col), val(val), lev(lev) {}

        friend bool operator<(const nonzero &a, const nonzero &b) {
            return a.col < b.col;
        }
    };

    struct sparse_vector {
        struct comp_indices {
            const std::deque<nonzero> &nz;

            comp_indices(const std::deque<nonzero> &nz) : nz(nz) {}

            bool operator()(int a, int b) const {
                return nz[a].col > nz[b].col;
            }
        };

        typedef std::priority_queue<int, std::vector<int>, comp_indices> priority_queue;

        int lfil;

        std::deque<nonzero>    nz;
        std::vector<ptrdiff_t> idx;
        priority_queue q;

        ptrdiff_t dia;

        sparse_vector(size_t n, int lfil)
            : lfil(lfil), idx(n, -1), q(comp_indices(nz)), dia(0) {}

        void add(ptrdiff_t col, const T &val, int lev) {
            if (idx[col] < 0) {
                if (lev <= lfil) {
                    int p = nz.size();
                    idx[col] = p;
                    nz.push_back(nonzero(col, val, lev));
                    if (col < dia) q.push(p);
                }
            } else {
                nonzero &a = nz[idx[col]];
                a.val += val;
                a.lev = std::min(a.lev, lev);
            }
        }

        typename std::deque<nonzero>::iterator begin() {
            return nz.begin();
        }

        typename std::deque<nonzero>::iterator end() {
            return nz.end();
        }

        nonzero& next_nonzero() {
            int p = q.top();
            q.pop();
            return nz[p];
        }

        void sort() {
            std::sort(nz.begin(), nz.end());
        }

        void reset(ptrdiff_t d) {
            for(const nonzero &e : nz) idx[e.col] = -1;
            nz.clear();
            dia = d;
        }
    };
};