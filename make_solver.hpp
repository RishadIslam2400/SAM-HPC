#pragma once

#include "CSRMatrix.hpp"

template <typename T, typename Precond, typename IterativeSolver>
class make_solver {
public:
    struct params {
        typename Precond::params precond;
        typename IterativeSolver::params solver;

        params() {}
    } prm;

    /** Sets up the preconditioner and creates the iterative solver. */
    make_solver(std::shared_ptr<CSRMatrix<T>> A, const params& prm = params())
        : prm(prm), n(A->m_rows), P(A, prm.precond), S(n, prm.solver) {}

    /** Computes the solution for the given system matrix \p A and the
     * right-hand side \p rhs.  Returns the number of iterations made and
     * the achieved residual as a ``std::tuple``. The solution vector
     * \p x provides initial approximation in input and holds the computed
     * solution on output.
     *
     * \rst
     * The system matrix may differ from the matrix used during
     * initialization. This may be used for the solution of non-stationary
     * problems with slowly changing coefficients. There is a strong chance
     * that a preconditioner built for a time step will act as a reasonably
     * good preconditioner for several subsequent time steps [DeSh12]_.
     * \endrst
     */
    std::tuple<size_t, T> operator()(const CSRMatrix<T>& A, const std::vector<T>& rhs, const std::vector<T>&& x) const {
        return S(A, P, rhs, x);
    }

    /** Computes the solution for the given right-hand side \p rhs.
     * Returns the number of iterations made and the achieved residual as a
     * ``std::tuple``. The solution vector \p x provides initial
     * approximation in input and holds the computed solution on output.
     */
    std::tuple<size_t, T> operator()(const std::vector<T> &rhs, std::vector<T> &&x) const {
        return S(P, rhs, x);
    }

    /** Acts as a preconditioner. That is, applies the solver to the
     * right-hand side \p rhs to get the solution \p x with zero initial
     * approximation.  Iterative methods usually use estimated residual for
     * exit condition.  For some problems the value of the estimated
     * residual can get too far from the true residual due to round-off
     * errors.  Nesting iterative solvers in this way may allow to shave
     * the last bits off the error. The method should not be used directly
     * but rather allows nesting ``make_solver`` classes as in the
     * following example:
     *
     * \rst
     * .. code-block:: cpp
     *
     *   typedef amgcl::make_solver<
     *     amgcl::make_solver<
     *       amgcl::amg<
     *         Backend, amgcl::coarsening::smoothed_aggregation, amgcl::relaxation::spai0
     *         >,
     *       amgcl::solver::cg<Backend>
     *       >,
     *     amgcl::solver::cg<Backend>
     *     > NestedSolver;
     * \endrst
     */
    void apply(const std::vector<T> &rhs, std::vector<T> &&x) const {
        x.clear();
        (*this)(rhs, x);
    }

    /// Returns const reference to the constructed preconditioner.
    const Precond& precond() const {
        return P;
    }

    /// Returns reference to the constructed preconditioner.
    Precond& precond() {
        return P;
    }

    /// Returns reference to the constructed iterative solver.
    const IterativeSolver& solver() const {
        return S;
    }

    /// Returns the system matrix in the backend format.
    std::shared_ptr<CSRMatrix<T>> system_matrix_ptr() const {
        return P.system_matrix_ptr();
    }

    CSRMatrix<T> const& system_matrix() const {
        return P.system_matrix();
    }

    /// Returns the size of the system matrix.
    size_t size() const {
        return n;
    }

    friend std::ostream& operator<<(std::ostream &os, const make_solver &p) {
        return os
            << "Solver\n======\n" << p.S << std::endl
            << "Preconditioner\n==============\n" << p.P;
    }

private:
    size_t n;
    Precond P;
    IterativeSolver S;
};