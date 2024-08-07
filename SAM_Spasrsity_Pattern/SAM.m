function [MM] = SAM(A_source,A_target, S)
% Computes sparse approximate map from a source matrix to a target matrix
% A_Source: Source matrix; typically any matrix in a sequence of linear
% system
% A_target: Target matrix; A matrix in the sequence of linear system that
% we want to map back to
% S: Boolean matrix; Given a priori sparsity pattern that we want to impose
% on the map

% Step 1: Preprocessing // Should we define another function
% precprocessing
% 1: Given sparsity pattern M and matrix A_source
% 2: maxSk = 0; maxRk = 0; {initialize max num of columns, max num of rows}
% 3: for k = 1 : n do {for each column do}
% 4:    sk = {i | (i, k) in M} {get indices; typically defined in advance}
% 5:    rk = null; {Initialize set of rows for kth LS problem}
% 6:    for all j in sk do
% 7:        t = find(aj) {find indices of nonzeros in column aj}
% 8:        rk = rk union t
% 9:    end for
% 10:   nnzk = #(sk) {#() gives number of elements in a set}
%       if nnzk > maxSk then maxSk = nnzk end if
%       if #(rk) > maxRk then maxRk = #(rk) end if
% 11: end for
% 12: Allocate maxRk × maxSk array for storing the LS matrices, maxRk vector for
% storing the right hand side, and maxSk vector for storing the solution.


% Step 2: Computing N
% 1: cnt = 0 {counts number of nonzeros in preconditioner}
% 2: (Preallocate space for A_tmp)
% 3: for k = 1 : n do
% 4:    A_tmp = A_source(rk, sk) {get submatrix indexed by rk and sk for LS problem}
% 5:    f = A_target (rk, k) {get rhs for LS problem}
% 6:    Solve LS A_tmp * n = f
% 7:    (possibly save residual, norm of residual, etc.)
% 8:    rowN[cnt + 1 : cnt + nnzk] = sk {assign indices in order of row ind. in sk}
% 9:    colN[cnt + 1 : cnt + nnzk] = k
% 10:   valN[cnt + 1 : cnt + nnzk] = n
% 11: end for
% 12: N = sparse(rowN; colN; valN) {convert into sparse matrix}

% Five parts of the function:
% 1. Finding indices I and J for each column
% 2. Calculating Local submatrix for each column from A_source
% 3. QR factorization of each local submatrix
% 4. Solving LS problem A^hat_j * n_j = a_0j
% 5. Assemble each column and output N

% 1. Preprocess the sparsity pattern
PP = logical(S); % Nonzero indices of each column
PP2 = logical(sparse(double(PP)*double(PP))); % Nonzero indices of each row for each entry in the column
nnzMM = nnz(PP2);
rowM = zeros(2*nnzMM,1);
colM = zeros(2*nnzMM,1);
valM = zeros(2*nnzMM,1);

% find the value of n
n = size(A_source, 2);

for j = 1:n
    % Initialize these later for efficiency...
    nz_M{j} = find(PP(:,j));
    nnz_M(j) = length(nz_M{j});
    nz_LS{j} = find(PP2(:,j));
    nnz_LS(j) = length(nz_LS{j});
end

max_col = max(nnz_M);
max_row = max(nnz_LS);
G = zeros(max_row,max_col);
M = zeros(max_row,1);
cntrM = 0;
for j = 1:n
    G(1:nnz_LS(j),1:nnz_M(j)) = A_source(nz_LS{j},nz_M{j});
    %[rC,cC] = size(G(1:nnz_LS(j),1:nnz_M(j)));
    M(1:nnz_M(j)) = G(1:nnz_LS(j),1:nnz_M(j))\A_target(nz_LS{j},j);
    rowM(cntrM+1:cntrM+nnz_M(j)) = nz_M{j};
    colM(cntrM+1:cntrM+nnz_M(j)) = j;
    valM(cntrM+1:cntrM+nnz_M(j)) = M(1:nnz_M(j));
    cntrM = cntrM+nnz_M(j);
end

MM = sparse(rowM(1:cntrM),colM(1:cntrM),valM(1:cntrM));

end

