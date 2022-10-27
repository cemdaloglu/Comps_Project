% Script containing all functions 
function [ x_reconstucted ] = l1_algorithm(A, b_arr)
    % returns reconstruction matrix of n x m, each column is one
    % reconstructed vector  for a sparsity from 1<=s<=m

    [m,n] = size(A);

    x_reconstucted = zeros(n, m);

    fprintf(1,'Sparsity Level: %3d\n', 0);
    for idx = 1:m
        fprintf(1,'\b\b\b%3.0f',idx);
        cvx_begin
            variable x_rec(n, 1)
            minimize(norm(x_rec, 1))
            subject to
                A * x_rec == b_arr(:, idx)
        cvx_end
        x_reconstucted(:, idx) = x_rec;
    end

    fprintf('\n');
   
end 







