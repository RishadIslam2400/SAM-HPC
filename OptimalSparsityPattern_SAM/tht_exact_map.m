% Load the data
load('data_early.mat');

% Create a directory to save sparsity patterns
sparsityDir = 'tht_sparsity_pattern';
if ~exist(sparsityDir, 'dir')
    mkdir(sparsityDir);
end

% Initialize a cell array to store the exact maps
thtMatrices = cell(1, numel(s));

for i = 1:numel(s)
    % Display the sparsity pattern of each sparse matrix
    A = K + s(i) * M;
    thtMatrices{i} = A;
    spy(A);
    title(['Sparsity Pattern of Matrix ', num2str(i)]);

    % Save the sparsity pattern as an image
    filename = fullfile(sparsityDir, ['sparsity_pattern_' num2str(i) '.png']);
    saveas(gcf, filename);
    
    % Close the current figure to avoid overlap
    close(gcf);
end

% Initialize a cell array to store the exact maps
exactMapArray = cell(1, numel(s) - 1);

% Create a directory to save sparsity patterns of inverse matrices
exactMapDir = 'tht_exact_map_sparsity_pattern';
if ~exist(exactMapDir, 'dir')
   mkdir(exactMapDir);
end

% Compute inverses and save sparsity patterns
for i = 1:numel(s)
    if i == 1
        A0 = thtMatrices{i};
    else
        A_k = thtMatrices{i};
        M_k = A_k\A0;
        exactMapArray{i-1} = M_k;

        % Show the sparsity pattern of the exact map
        spy(M_k);
        title(['Exact Map (step' num2str(i) '.mat)']);
        xlabel('Column Index');
        ylabel('Row Index');

        % Save the sparsity pattern plot of the exact map
        exactMapFilename = fullfile(exactMapDir, [num2str(i) '.png']);
        saveas(gcf, exactMapFilename);
        close(gcf);
    end
end

% Save inverseMatArray to a .mat file
save('tht_matrices.mat', 'thtMatrices', '-v7.3', '-nocompression');
save('tht_exact_maps.mat', 'exactMapArray', '-v7.3', '-nocompression');
