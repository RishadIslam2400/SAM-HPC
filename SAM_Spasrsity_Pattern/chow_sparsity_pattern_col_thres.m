function [S] = chow_sparsity_pattern_col_thres(A, tau)
% Returns the boolean matrix representing the sparsity pattern of the input
% matrix according to method given by Edmond Chow
%   

% Get the size of A
[m, n] = size(A);

% Create the sparsity pattern
S = sparse(m, n);

for j = 1:n
    % Compute the threshold for the current column
    thresh = (1 - tau) * max(abs(A(:, j)));

    % Find the row indices where the absolute value is greater than the threshold
    idx = abs(A(:, j)) > thresh;

    % Set the corresponding entries in the sparsity pattern matrix S to 1
    S(idx, j) = 1;
end

% Set the diagonal entries to one
S(1:size(A,1)+1:end) = 1;

