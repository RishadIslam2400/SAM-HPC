% Load inverse matrices
exact_maps = load('bubble_exact_maps.mat');
exact_maps = exact_maps.exactMapArray;

% Create a directory to save sparse exact maps
exactmap_dir = 'bubble_exact_map_lfil_7';
if ~exist(exactmap_dir, 'dir')
   mkdir(exactmap_dir);
end

for i = 1:numel(exact_maps)
    A_k = exact_maps{i};
    S_k = chow_sparsity_pattern_lfil(A_k, 7);

    % Show the sparsity pattern of matrix A
    spy(S_k);
    title(['Sparse Exact map Step' num2str(i+1)]);
    % xlabel('Column Index');
    % ylabel('Row Index');

    % Save the sparsity pattern plot
    sparsityFilename = fullfile(exactmap_dir, ['step' num2str(i+1) '.png']);
    saveas(gcf, sparsityFilename);
    close(gcf);
end