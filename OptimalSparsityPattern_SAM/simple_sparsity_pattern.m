function [S] = simple_sparsity_pattern(A)
% Returns the boolean matrix representing the sparsity pattern of the input
% matrix
%   
[I, J] = find(A);
S = sparse(I, J, 1);

% Level 2 neighbors
% S = S * S;
end

