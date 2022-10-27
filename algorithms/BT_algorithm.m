function [ x_reconstructed ] = BT_algorithm(A, b_arr)
    [m,n] = size(A);
    x_reconstructed = zeros(n, m);
    
    for idx = 1:m
        cor = abs(A' * b_arr(:, idx));
        [a, b] = maxk(cor, idx);
    
        cvx_begin
            variable x_rec(idx, 1)
            minimize(norm(b_arr(:, idx) - A(:, b) * x_rec, 2))
        cvx_end
    
        x_sup = zeros(n, 1);
        x_sup(b, :) = x_rec;
    
        x_reconstructed(:, idx) = x_sup;
    end
end 