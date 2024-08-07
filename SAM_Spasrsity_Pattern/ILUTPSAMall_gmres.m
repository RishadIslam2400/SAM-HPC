%Computes an ILUTP preconditioner for the first system in the sequence and
%updates for all others using the SAM update

% Parameters: 
% nx: number of grids along x direction
% ny: number of grids along y direction
% n: number of column/rows in the discretized matrix (nx * ny)
% u0: inital value (arbirtray matrix)

function [iterations, soltime, prectime,fdtime,gmresinfo,back,sol] = ILUTPSAMall_gmres(nx,ny,n,u0)

dx = 1/(nx+1); % Step size along x direction
dy = 1/(ny+1); % Step size along y direction
x0 = 0; % initial boundary condition
y0 = 0; % initial boundary condition


stopIter = 1000; % Max iteration

% Containers for saving result of each iteration
% Return values of the function
iterations = zeros(stopIter,1);
soltime = zeros(stopIter,1);
prectime = zeros(stopIter,1);
fdtime = cell(stopIter,2);
gmresinfo = cell(stopIter,2);
sol = cell(stopIter,1);
back = 0;


u = u0;


ftol = 1.e-4; 
p = inf;
k = 0; % iteration count
F = inf;
alpha = 1.e-4;
max_m = 10;
rtol = 5.e-10;
atol = 5.e-10;
droptol = 1.e-5; lfil = 10;
xx0 = zeros(n,1); ttol = 1.e-10;  max_it = n; restart = 50;

% Stopping condition for the sequence of linear system
% Means we achieved a good enough  approximation
while  (norm(F) > (ftol*rtol + atol)) && (k <= stopIter)
    k = k + 1
    norm(F)
    tic
    % Calculate the sequence of parametrized linear systems from nonlienear
    % PDEs
    [F,Jac,~,~] = cd2d_nonlinear(nx,ny,dx,dy,x0,y0,u, ...
        @pcoef,@qcoef,@pcoefdx,@qcoefdx,@rcoef,@scoef,@tcoef,@fcoef, ...
        @sbc,@wbc,@nbc,@ebc);
    fdtime{k,1} = toc;
    if k == 1
        ftol = norm(F);
    end

    % Sparsity pattern of matrix A/Jac
%     figure;
%     spy(Jac);
%     filename = sprintf('images/sparsity_pattern_%d.png', k);
%     saveas(gcf,filename);
%     close(gcf);

    % Save the Jac matrix
%     mat_filename = sprintf('data/matrix_%d.mat', k);
%     save(mat_filename, 'Jac');

    % During the first iteration compute the ILUTP preconditioner for the
    % first linear system
    % Also compute the sparsity pattern of the first linear system    
    if k == 1
        tic
        [L,U,perm,~,Jp] = my_ILUTP(Jac,droptol,lfil);
        prectime(k) = toc;
        mats = {'ILU',Jp,L,U};
        %[I,J] = find(Jac); % Anything i have to do with the sparsity pattern has to be done here insidethe find, PP1 and PP2 just takes columns and rows
        %findJ = sparse(I,J,1);
        %findJ = simple_sparsity_pattern(Jac);
        %findJ = chow_sparsity_pattern_global_thres(Jac, 0.1);
%         PP = logical(findJ); % Nonzero indices of each column
%         PP2 = logical(sparse(double(PP)*double(PP))); % Nonzero indices of each row for each entry in the column
%         nnzMM = nnz(PP2);
%         rowM = zeros(2*nnzMM,1);
%         colM = zeros(2*nnzMM,1);
%         valM = zeros(2*nnzMM,1);
        J0 = Jac; % Taking the first matrix in the sequence as A_0
    else
        % Compute SAM
        % During the second iteration: Compute SAM using the sparsity
        % pattern computed in the first iteration
        tic
        if k == 2 % Preprocess the sparsity pattern (here, we use patt(A))
%             for j = 1:n
%             %Initialize these later for efficiency...
%                 nz_M{j} = find(PP(:,j));
%                 nnz_M(j) = length(nz_M{j});
%                 nz_LS{j} = find(PP2(:,j));
%                 nnz_LS(j) = length(nz_LS{j});
%             end
%             
%             max_col = max(nnz_M);
%             max_row = max(nnz_LS);
%             G = zeros(max_row,max_col);
%             M = zeros(max_row,1);
%             cntrM = 0;
%             for j = 1:n
%                G(1:nnz_LS(j),1:nnz_M(j)) = Jac(nz_LS{j},nz_M{j});
%                [rC,cC] = size(G(1:nnz_LS(j),1:nnz_M(j)));
%                M(1:nnz_M(j)) = G(1:nnz_LS(j),1:nnz_M(j))\J0(nz_LS{j},j);
%                rowM(cntrM+1:cntrM+nnz_M(j)) = nz_M{j};
%                colM(cntrM+1:cntrM+nnz_M(j)) = j;
%                valM(cntrM+1:cntrM+nnz_M(j)) = M(1:nnz_M(j));
%                cntrM = cntrM+nnz_M(j);
%             end
%             MM = sparse(rowM(1:cntrM),colM(1:cntrM),valM(1:cntrM));
            findJ = simple_sparsity_pattern(Jac);
            %findJ = chow_sparsity_pattern_global_thres(Jac, 0.1);
            %findJ = chow_sparsity_pattern_col_thres(Jac, 0.5);
            %findJ = chow_sparsity_pattern_lfil(Jac, 3);
            MM = SAM(Jac, J0, findJ);
            % permute MM (with perm)
            %                 MM = MM(:,perm(:));
            
            % Here we directly permute the jacobian; you can instead permute
            % the map (see immediately above) OR if you use ILU/ILUT you
            % can remove any use of the permutation array
            mats = {'SAM',Jac(:,perm(:)),L,U,MM};
            % In other iterations the preprocessing for allocating space is
            % not done
        elseif k > 2
            
%             max_col = max(nnz_M);
%             max_row = max(nnz_LS);
%             G = zeros(max_row,max_col);
%             M = zeros(max_row,1);
%             cntrM = 0;
%             for j = 1:n
%                G(1:nnz_LS(j),1:nnz_M(j)) = Jac(nz_LS{j},nz_M{j});
%                %[rC,cC] = size(G(1:nnz_LS(k),1:nnz_M(k)));
%                M(1:nnz_M(j)) = G(1:nnz_LS(j),1:nnz_M(j))\J0(nz_LS{j},j);
%                rowM(cntrM+1:cntrM+nnz_M(j)) = nz_M{j};
%                colM(cntrM+1:cntrM+nnz_M(j)) = j;
%                valM(cntrM+1:cntrM+nnz_M(j)) = M(1:nnz_M(j));
%                cntrM = cntrM+nnz_M(j);
%             end
%             MM = sparse(rowM(1:cntrM),colM(1:cntrM),valM(1:cntrM));
            findJ = simple_sparsity_pattern(Jac);
            %findJ = chow_sparsity_pattern_global_thres(Jac, 0.1);
            %findJ = chow_sparsity_pattern_col_thres(Jac, 0.5);
            %findJ = chow_sparsity_pattern_lfil(Jac, 3);
            MM = SAM(Jac, J0, findJ);
            % permute MM (with perm)
            MM = MM(:,perm(:));
            
            mats = {'SAM',Jac,L,U,MM};
        end
        prectime(k) = toc;
    end
    
    tic
    [vtemp,r,r_nrm,iter,~] = mgmres(mats,F,xx0,ttol,max_it,restart,'Right');
    if k == 1
        p = (U\(L\vtemp));
    elseif k > 1
        p = MM*(U\(L\vtemp));
    end
    soltime(k) = toc;
    iterations(k) = iter;
    gmresinfo{k,1} = r;
    gmresinfo{k,2} = r_nrm;
    tic
    [p,back] = backstep(back,u,nx,ny,dx,dy,x0,y0,max_m,F,p,alpha);
    fdtime{k,2} = toc;
    u = u + p;
    sol{k} = u;
end