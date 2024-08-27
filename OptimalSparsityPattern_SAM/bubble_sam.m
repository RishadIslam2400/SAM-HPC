% load the data
data = load("sequence_of_A.mat");
sequence_of_A = data.matArray;

% Initialize a cell array to store the exact maps
approxMapArray = cell(1, 20);

% Create a directory to save sparsity patterns of inverse matrices
approxmap_dir = 'bubble_approx_map_global_sparsity_pattern_squared_0.0001';
if ~exist(approxmap_dir, 'dir')
   mkdir(approxmap_dir);
end

% Compute inverses and save sparsity patterns
for i = 1:21
    if i == 1
        A0 = sequence_of_A{i};
    else
        A_k = sequence_of_A{i};
        S_k = chow_sparsity_pattern_global_thres(A_k, 0.0001);
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

% Save inverseMatArray to a .mat file
save('bubble_approx_maps_global_squared_0.0001.mat', 'approxMapArray');

% % Initialize a cell array to store the exact maps
% approxMapArray_global = cell(1, 20);
% 
% % Create a directory to save sparsity patterns of inverse matrices
% approxmap_global = 'approx_map_global_sparsity_pattern';
% if ~exist(approxmap_global, 'dir')
%    mkdir(approxmap_global);
% end
% 
% % Compute inverses and save sparsity patterns
% for i = 1:21
%     if i == 1
%         A0 = sequence_of_A{i};
%     else
%         A_k = sequence_of_A{i};
%         S_k = chow_sparsity_pattern_global_thres(A_k, 0.1);
%         M_k = SAM(A_k, A0, S_k);
% 
%         approxMapArray_global{i-1} = M_k;
% 
%         % Show the sparsity pattern of the approximate map
%         spy(M_k);
%         title(['Approximate Map (step' num2str(i) '.mat)']);
%         xlabel('Column Index');
%         ylabel('Row Index');
% 
%         % Save the sparsity pattern plot of the exact map
%         approxMapFilename = fullfile(approxmap_global, ['step' num2str(i) '.png']);
%         saveas(gcf, approxMapFilename);
%         close(gcf);
%     end
% end
% 
% % Save inverseMatArray to a .mat file
% save('approx_maps_global_sparsity_pattern.mat', 'approxMapArray_global');
% 
% % Initialize a cell array to store the exact maps
% approxMapArray_col = cell(1, 20);
% 
% % Create a directory to save sparsity patterns of inverse matrices
% approxmap_col = 'approx_map_col_sparsity_pattern';
% if ~exist(approxmap_col, 'dir')
%    mkdir(approxmap_col);
% end
% 
% % Compute inverses and save sparsity patterns
% for i = 1:21
%     if i == 1
%         A0 = sequence_of_A{i};
%     else
%         A_k = sequence_of_A{i};
%         S_k = chow_sparsity_pattern_col_thres(A_k, 0.5);
%         M_k = SAM(A_k, A0, S_k);
% 
%         approxMapArray_col{i-1} = M_k;
% 
%         % Show the sparsity pattern of the approximate map
%         spy(M_k);
%         title(['Approximate Map (step' num2str(i) '.mat)']);
%         xlabel('Column Index');
%         ylabel('Row Index');
% 
%         % Save the sparsity pattern plot of the exact map
%         approxMapFilename = fullfile(approxmap_col, ['step' num2str(i) '.png']);
%         saveas(gcf, approxMapFilename);
%         close(gcf);
%     end
% end
% 
% % Save inverseMatArray to a .mat file
% save('approx_maps_col_sparsity_pattern.mat', 'approxMapArray_col');
% 
% % Initialize a cell array to store the exact maps
% approxMapArray_lfil = cell(1, 20);
% 
% % Create a directory to save sparsity patterns of inverse matrices
% approxmap_lfil = 'approx_map_lfil_sparsity_pattern';
% if ~exist(approxmap_lfil, 'dir')
%    mkdir(approxmap_lfil);
% end
% 
% % Compute inverses and save sparsity patterns
% for i = 1:21
%     if i == 1
%         A0 = sequence_of_A{i};
%     else
%         A_k = sequence_of_A{i};
%         S_k = chow_sparsity_pattern_lfil(A_k, 5);
%         M_k = SAM(A_k, A0, S_k);
% 
%         approxMapArray_lfil{i-1} = M_k;
% 
%         % Show the sparsity pattern of the approximate map
%         spy(M_k);
%         title(['Approximate Map (step' num2str(i) '.mat)']);
%         xlabel('Column Index');
%         ylabel('Row Index');
% 
%         % Save the sparsity pattern plot of the exact map
%         approxMapFilename = fullfile(approxmap_lfil, ['step' num2str(i) '.png']);
%         saveas(gcf, approxMapFilename);
%         close(gcf);
%     end
% end
% 
% % Save inverseMatArray to a .mat file
% save('approx_maps_lfil_sparsity_pattern.mat', 'approxMapArray_lfil');