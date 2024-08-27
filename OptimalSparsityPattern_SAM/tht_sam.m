load('tht_matrices.mat');

approxMapArray = cell(1, numel(thtMatrices) - 1);
approxmap_dir = 'tht_approx_map_simple_sparsity_pattern';
if ~exist("approxmap_dir", 'dir')
   mkdir(approxmap_dir);
end

for i = 1:numel(thtMatrices)
    if i == 1
        A0 = thtMatrices{i};
    else
        A_k = thtMatrices{i};
        S_k = simple_sparsity_pattern(A_k);
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

save('tht_approx_maps_simple.mat', 'approxMapArray');