load('cd2d_matrices.mat');

approxMapArray = cell(1, numel(matrices) - 1);
approxmap_dir = 'cd2d_approx_maps_column_sparsity_pattern_squared_0.8';
if ~exist("approxmap_dir", 'dir')
   mkdir(approxmap_dir);
end

for i = 1:numel(matrices)
    if i == 1
        A0 = matrices{i};
    else
        A_k = matrices{i};
        S_k = chow_sparsity_pattern_col_thres(A_k, 0.8);
        M_k = SAM(A_k, A0, S_k);

        approxMapArray{i-1} = M_k;

        % Show the sparsity pattern of the approximate map
        spy(M_k);
        title(['Approximate Map (step' num2str(i) '.mat)']);
        % xlabel('Column Index');
        % ylabel('Row Index');

        % Save the sparsity pattern plot of the exact map
        approxMapFilename = fullfile(approxmap_dir, ['step' num2str(i) '.png']);
        saveas(gcf, approxMapFilename);
        close(gcf);
    end
end

save('cd2d_approx_maps_column_squared_0.8.mat', 'approxMapArray');