function [S] = chow_sparsity_pattern_lfil(A, lfil)
% Returns the boolean matrix representing the sparsity pattern of the input
% matrix according to method given by Edmond Chow
%   

% Get the size of A
[m, n] = size(A);

% Initialize row, column, and value arrays to construct the sparse matrix S
max_elements = lfil * n; % maximum number of elements in S
row_idx = zeros(max_elements, 1);
col_idx = zeros(max_elements, 1);
values = ones(max_elements, 1);

% Initialize index variable for tracking current position in arrays
idx = 1;

for j = 1:n
    % Find the row indices of non-zero entries in the current column
    idx_col = find(A(:, j));
    
    % If lfil is greater than the number of non-zero entries, take all non-zero entries
    if lfil >= numel(idx_col)
        lfil_col = idx_col;
    else
        % Take the lfil largest non-zero elements
        [~, sorted_idx] = sort(abs(A(idx_col, j)), 'descend');
        lfil_col = idx_col(sorted_idx(1:lfil));
    end
    
    % Update row_idx and col_idx arrays with the selected entries
    num_selected = numel(lfil_col);
    row_idx(idx:idx+num_selected-1) = lfil_col;
    col_idx(idx:idx+num_selected-1) = j;
    
    % Update index variable
    idx = idx + num_selected;
end

% Construct the sparsity pattern matrix S
S = sparse(row_idx(1:idx-1), col_idx(1:idx-1), values(1:idx-1), m, n);

% Set the diagonal entries to one
S(1:size(A,1)+1:end) = 1;

S = S * S;

