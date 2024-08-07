function [S] = chow_sparsity_pattern_global_thres(A, thresh)
% Returns the boolean matrix representing the sparsity pattern of the input
% matrix according to method given by Edmond Chow
%   
% Input:    A: Input matrix to create the sparsity pattern
%           thresh: Global threshold
% Output:   S: Sparsity pattern according to the formula
%           A_ij = 1 if i = j or |D^(-1/2)AD^(-1/2)| > thresh
%                = 0 otherwise
%           D_ii = |A_ii| if |A_ii| > 0
%                = 1 otherwise

% Set up the diagonal matrix
diagonal_values = diag(abs(A));
diagonal_values(diagonal_values == 0) = 1;
%D = spdiags(abs(diagonal_values), 0, size(A, 1), size(A, 2));

% Normalize the matrix A by the diagonal matrix
D_inv = spdiags(1./sqrt(diagonal_values), 0, size(A, 1), size(A, 2)); 
A_norm = D_inv * A * D_inv;

% Create the sparsity pattern
[I, J, V] = find(A_norm);
S = sparse(I, J, V, size(A, 1), size(A, 2));
S(1:size(A, 1) + 1: end) = 1;
S(abs(S) > thresh) = 1;
end


