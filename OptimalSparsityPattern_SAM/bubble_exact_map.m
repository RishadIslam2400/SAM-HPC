% Load the data
directory = 'data/';

% Initialize the matArray structure
matArray = cell(1, 21);

% Create a directory to save sparsity patterns
sparsityDir = 'bubble_sparsity_pattern';
if ~exist(sparsityDir, 'dir')
    mkdir(sparsityDir);
end

% Iterate through each .mat file and load its contents
for i = 1:21
    filename = ['step' num2str(i) '.mat'];
    filePath = fullfile(directory, filename);
    
    % Load the .mat file
    matData = load(filePath);
    
    % Extract coefficient matrix A and right-hand side vector b
    A = matData.A;
    
    % Store the data in the matArray structure
    matArray{i} = A;

    % Show the sparsity pattern of matrix A
    spy(A);

    % Customize the plot
    ax = gca; % Get current axes
    ax.FontSize = 20; % Set font size
    ax.LineWidth = 2; % Set box line width

    % Save the sparsity pattern plot
    sparsityFilename = fullfile(sparsityDir, ['step' num2str(i) '.png']);
    saveas(gcf, sparsityFilename);
    
    % Save the sparsity pattern plot as FIG
    sparsityFilenameFIG = fullfile(sparsityDir, ['step' num2str(i) '.fig']);
    savefig(sparsityFilenameFIG);

    close(gcf);
end

% Save matArray to a .mat file
% save('sequence_of_A.mat', 'matArray');

% Create a directory to save sparsity patterns of inverse matrices
% exactMapDir = 'exact_map_sparsity_pattern';
% if ~exist(exactMapDir, 'dir')
%    mkdir(exactMapDir);
% end

% Initialize a cell array to store the exact maps
% exactMapArray = cell(1, 20);

% Compute inverses and save sparsity patterns
% for i = 1:21
%     if i == 1
%         A0 = matArray{i};
%     else
%         A_k = matArray{i};
%         M_k = A_k\A0;
%         exactMapArray{i-1} = M_k;
% 
%         % Show the sparsity pattern of the exact map
%         % spy(M_k);
%         % title(['Exact Map (step' num2str(i) '.mat)']);
%         % xlabel('Column Index');
%         % ylabel('Row Index');
% 
%         % Save the sparsity pattern plot of the exact map
%         % exactMapFilename = fullfile(exactMapDir, [num2str(i) '.png']);
%         % saveas(gcf, exactMapFilename);
%         % close(gcf);
%     end
% end

% Save inverseMatArray to a .mat file
% save('exact_maps.mat', 'exactMapArray', '-v7.3', '-nocompression');
