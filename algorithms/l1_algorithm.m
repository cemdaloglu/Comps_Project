% Script containing all functions 
function [ x_reconstucted ] = l1_algorithm(A, b_arr)
    % returns reconstruction matrix of n x m, each column is one
    % reconstructed vector  for a sparsity from 1<=s<=m

    [m,n] = size(A);

    x_reconstucted = zeros(n, m);
    b_arr = b_arr';
    
    fprintf(1,'Sparsity Level: %3d\n', 0);
    for idx = 1:m
        fprintf(1,'\b\b\b%3.0f',idx);
        
        cvx_begin 
            variable x_rec(n, 1)
            minimize(norm(x_rec, 1))
            subject to
                A * x_rec == b_arr(:, idx);
        cvx_end
        x_reconstucted(:, idx) = x_rec;
    end
    fprintf('\n');
   
end 


function [ x_rec_random, x_rec_fourier ] = OMP_algorithm(A_random, b_random, A_fourier, b_fourier)
end 


function [ x_rec_random, x_rec_fourier ] = MP_algorithm(A_random, b_random, A_fourier, b_fourier)
    x_rec_random_0 = zeros(size(A_random,2), 1);    % initialize 
    x_rec_fourier_0 = zeros(size(A_fourier,2), 1);
    r_rand_0 = b_random; 
    r_fourier_0 = b_fourier;
    k = 0; 
    scalarproducts = A_random'*x;
    residual = x-scalarproducts(3)*dictionary(:,3);
    scalarproducts = dictionary(:,1:2)'*residual;
end 


function [ x_rec_random, x_rec_fourier ] = IHT_algorithm(A_random, b_random, A_fourier, b_fourier, sparsity_level)
end 


function [ x_rec_random, x_rec_fourier ] = CoSaMP_algorithm(A_random, b_random, A_fourier, b_fourier, sparsity_level)
end 


function [ x_rec_random, x_rec_fourier ] = BT_algorithm(A_random, b_random, A_fourier, b_fourier, sparsity_level)
end 


function [ x_rec_random, x_rec_fourier ] = HTP_algorithm(A_random, b_random, A_fourier, b_fourier)
end 


function [ x_rec_random, x_rec_fourier ] = SP_algorithm(A_random, b_random, A_fourier, b_fourier, sparsity_level)
end 

