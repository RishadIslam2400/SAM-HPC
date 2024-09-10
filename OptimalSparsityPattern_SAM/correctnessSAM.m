% load data
directory = 'data/';

% Initialize the matArray structure
matArray = cell(1, 70);


% Iterate through each .mat file and load its contents
for i = 1:70
    filename = ['matrix_' num2str(i) '.mat'];
    filePath = fullfile(directory, filename);
    
    % Load the .mat file
    matData = load(filePath);
    
    % Extract coefficient matrix A and right-hand side vector b
    A = matData.Jac;
    
    % Store the data in the matArray structure
    matArray{i} = A;
end

% Run the SAM algorithm
for i = 1:2 % change it to 70 for running all the matrices
    if i == 1
        A0 = matArray{i};
    else
        A_k = matArray{i};
        S_k = simple_sparsity_pattern(A0);
        M_k = SAM(A_k, A0, S_k);

       % Not doing anything with M_k yet
    end
end