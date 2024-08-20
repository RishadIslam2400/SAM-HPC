import numpy as np
from scipy.sparse import csc_matrix, save_npz, load_npz
from scipy.linalg import qr
# # Matrix 1: py_matrix_10.csc
# values_1 = [0.01, 1, 0.1, 1.1, 2.1, 1.2, 2.2, 3.2, 2.3, 3.3, 4.3, 3.4, 4.4, 5.4, 4.5, 5.5, 6.5, 5.6, 6.6, 7.6, 6.7, 7.7, 8.7, 7.8, 8.8, 9.8, 8.9, 9.9]
# row_indices_1 = [0, 1, 0, 1, 2, 1, 2, 3, 2, 3, 4, 3, 4, 5, 4, 5, 6, 5, 6, 7, 6, 7, 8, 7, 8, 9, 8, 9]
# col_pointers_1 = [0, 2, 5, 8, 11, 14, 17, 20, 23, 26, 28]

# # Create the CSC matrix for the first one
# matrix_1 = csc_matrix((values_1, row_indices_1, col_pointers_1), shape=(10, 10))

# # Save the first matrix as a .npz file
# save_npz("py_matrix_10.npz", matrix_1)

# # Matrix 2: py_matrixA_10.csc
# values_2 = [0.24, 0.67, 0.18, 1.23, 0.56, 0.89, 1.34, 0.75, 0.98, 0.12, 0.45, 0.36, 0.83, 0.59, 0.26]
# row_indices_2 = [1, 0, 5, 2, 3, 6, 7, 9, 0, 4, 5, 2, 3, 7, 8]
# col_pointers_2 = [0, 1, 3, 5, 6, 8, 10, 11, 12, 13, 15]

# # Create the CSC matrix for the second one
# matrix_2 = csc_matrix((values_2, row_indices_2, col_pointers_2), shape=(10, 10))

# # Save the second matrix as a .npz file
# save_npz("py_matrixA_10.npz", matrix_2)

# print("Matrices saved as .npz files.")

# values = [0.8364, 0.4329, 0.9683, 0.3754, 0.8913, 0.2937, 0.7618, 0.5129, 0.6403, 0.3411, 
#           0.8922, 0.1104, 0.2755, 0.5893, 0.3652, 0.7562, 0.9132, 0.1223, 0.8346, 0.7631, 
#           0.2839, 0.4984, 0.4321, 0.6734, 0.7541, 0.1239, 0.9834, 0.3645, 0.7451, 0.8321]

# row_indices = [0, 7, 1, 5, 4, 3, 6, 2, 9, 8, 
#                0, 6, 2, 3, 4, 1, 8, 7, 2, 5, 
#                6, 1, 0, 9, 3, 7, 8, 4, 5, 0]

# col_pointers = [0, 3, 6, 8, 11, 15, 17, 21, 24, 27, 30]

# # Create the CSC matrix
# csc_matrix_data = csc_matrix((values, row_indices, col_pointers), shape=(10, 10))

# # Save the matrix to a .npz file
# save_npz("py_matrixb_10.npz", csc_matrix_data)





def find_max_length(my_dict):
    max_length = 0 
    for key, value in my_dict.items():
        if len(value) > max_length:
            max_length = len(value)
    
    return max_length

def remove_duplicates(my_dict):
    for key, value in my_dict.items():
        seen = set()  
        my_dict[key] = [x for x in value if not (x in seen or seen.add(x))]
    
    return my_dict


def find_actual_size(my_dict):
    size_dict = {}
    for key, value in my_dict.items():
        size_dict[key] = len(value)
    
    return size_dict

N = 10
# Load the first matrix
mat_S = load_npz("py_matrix_10.npz")
mat_S = mat_S.toarray()

# Load the second matrix
mat_A = load_npz("py_matrixA_10.npz")
mat_A = mat_A.toarray()

mat_b = load_npz("py_matrixb_10.npz")
mat_b = mat_b.toarray()

mat_S = (mat_S != 0).astype(int)
print(mat_S)




Sk_dict = {}
Rk_dict = {}
submat_dict = {}
for k in range(N):
    Sk_dict[k] = []
    Rk_dict[k] = []
    


for k in range(N):
    for row, element in enumerate(mat_S[:, k]):
        if element == 1:
            Sk_dict[k].append(row)

print(Sk_dict)
max_Sk = find_max_length(Sk_dict)

print(mat_A)
for k, indices in Sk_dict.items():
    for sub_mat_col_idx in indices:
        for row, element in enumerate(mat_A[:, sub_mat_col_idx]):
            if mat_A[row,sub_mat_col_idx] != 0.0:
                Rk_dict[k].append(row)


Rk_dict = remove_duplicates(Rk_dict)
max_Rk = find_max_length(Rk_dict)
print(Rk_dict)

Rk_actual_sizes = find_actual_size(Rk_dict)
Sk_actual_sizes = find_actual_size(Sk_dict)

print(f"Max column number in the submatrix: {max_Sk}")
print(f"Max row number in the submatrix: {max_Rk}")
print(Sk_actual_sizes)
print(Rk_actual_sizes)

################################################################################
sub_mat_dict = {}
b_dict = {}

print(mat_b)

for k in range(N):
    row_idx_k = Rk_dict[k]
    col_idx_k = Sk_dict[k]
    
    submatrix = np.zeros((max_Rk, max_Sk))
    

    for i, row_idx in enumerate(row_idx_k):
        for j, col_idx in enumerate(col_idx_k):
            if row_idx < mat_A.shape[0] and col_idx < mat_A.shape[1]:  
                submatrix[i, j] = mat_A[row_idx, col_idx]

    # Save the submatrix in the dictionary
    sub_mat_dict[k] = submatrix
    b_dict[k] = []

    for row_idx in row_idx_k:
        b_dict[k].append(float(mat_b[row_idx, k]))

# Print the submatrices to verify
for k, submatrix in sub_mat_dict.items():
    print(f"Submatrix for k={k}:")
    print(submatrix)

print(b_dict)

    
# Store the Q and R matrices for each submatrix
Q_dict = {}
R_dict = {}

for k, submatrix in sub_mat_dict.items():
    # Get the actual row and column sizes for this submatrix
    actual_rows = Rk_actual_sizes[k]  # Actual number of rows
    actual_cols = Sk_actual_sizes[k]  # Actual number of columns

    # Slice the submatrix to only include the actual rows and columns
    actual_submatrix = submatrix[:actual_rows, :actual_cols]

    # Perform QR decomposition on the actual submatrix
    Q, R = np.linalg.qr(actual_submatrix)  # Perform QR decomposition
    
    # Store the Q and R matrices in dictionaries
    Q_dict[k] = Q
    R_dict[k] = R
    
    # Print the Q and R matrices to verify
    print(f"Q matrix for submatrix k={k}:")
    print(Q)
    print(f"R matrix for submatrix k={k}:")
    print(R)