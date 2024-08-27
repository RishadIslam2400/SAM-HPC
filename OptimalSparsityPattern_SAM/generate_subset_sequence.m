% Define the size of the matrix
matrix_size = 1080;

% Number of matrices in the sequence
num_matrices = 11;

% Generate sequence of matrices 
subset_matrices = cell(1, num_matrices);

sparsity_level = 0.0001;

% Generate the first matrix with random sparsity
matrices{1} = sprandn(matrix_size, matrix_size, sparsity_level);

% Ensure non-zero diagonal elements
matrices{1}(1:matrix_size+1:end) = rand(matrix_size, 1) + 0.1;

% Generate subsequent matrices with increasing sparsity
for i = 2:num_matrices
    % Sparsity increases gradually
    sparsity_level = sparsity_level + 0.00001;
    
    % Ensure subset order by copying the previous matrix and increasing sparsity
    matrices{i} = matrices{i-1};
    
    % Introduce additional random sparsity
    additional_sparsity = sprandn(matrix_size, matrix_size, sparsity_level);
    
    % Add the additional sparsity to the current matrix
    matrices{i} = max(matrices{i}, additional_sparsity);
    
    % Ensure non-zero diagonal elements
    matrices{i}(1:matrix_size+1:end) = rand(matrix_size, 1) + 0.1;
end

% Create a folder to save sparsity patterns
sparsity_folder = 'subset_sparsity_patterns';
if ~exist(sparsity_folder, 'dir')
    mkdir(sparsity_folder);
end

% Display the sparsity pattern of each matrix and save the image
for i = 1:num_matrices
    % Display sparsity pattern using spy function
    figure;
    spy(matrices{i});
    title(sprintf('Sparsity Pattern - Matrix %d', i));
    
    % Save the image in the folder
    sparsity_img_filename = fullfile(sparsity_folder, sprintf('subset_sparsity_pattern_%d.png', i));
    saveas(gcf, sparsity_img_filename);
    close(gcf); % Close the figure to prevent accumulation
end

% Save the sequence of matrices and corresponding right-hand sides in a .mat file
sequence_filename = 'subset_matrix_sequence.mat';
save(sequence_filename, 'matrices');