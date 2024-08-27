% load the data
data = load("subset_matrix_sequence.mat");
subset_sequence_of_A = data.matrices;


approxMapArray = cell(1, numel(subset_sequence_of_A) - 1);
approxmap_dir = 'subset_approx_maps_col_sparsity_pattern';
if ~exist("approxmap_dir", 'dir')
   mkdir(approxmap_dir);
end

for i = 1:numel(subset_sequence_of_A)
    if i == 1
        A0 = subset_sequence_of_A{i};
    else
        A_k = subset_sequence_of_A{i};
        S_k = chow_sparsity_pattern_col_thres(A_k, 0.5);
        M_k = SAM(A_k, A0, S_k);

        approxMapArray{i-1} = M_k;

        % Show the sparsity pattern of the approximate map
        spy(M_k);
        title(['Approximate Map (step' num2str(i) '.mat)']);
        xlabel('Column Index');
        ylabel('Row Index');

        % Save the sparsity pattern plot of the exact map
        approxMapFilename = fullfile(approxmap_dir, ['step' num2str(i) '.png']);
        saveas(gcf, approxMapFilename);
        close(gcf);
    end
end

save('subset_approx_maps_col.mat', 'approxMapArray');