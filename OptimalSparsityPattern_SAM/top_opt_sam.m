load('top_opt_matrices.mat');

approxMapArray = cell(1, numel(Mats) - 1);
approxmap_dir = 'top_opt_approx_maps_init_sparsity_pattern';
if ~exist("approxmap_dir", 'dir')
   mkdir(approxmap_dir);
end

for i = 1:numel(Mats)
    if i == 1
        A0 = Mats{i};
    else
        A_k = Mats{i};
        S_k = simple_sparsity_pattern(A0);
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

save('top_opt_approx_maps_init.mat', 'approxMapArray');