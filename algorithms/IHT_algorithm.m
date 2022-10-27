function [ x_reconstructed ] = IHT_algorithm(A, b_arr)

    [m,n] = size(A);
    x_reconstructed = zeros(n, m);
    
    for idx = 1:m
        res = b_arr(:, idx);
        count = 0;
        x_rec = zeros(n, 1);
        while norm(res, 2) > 10^(-6) && count < m
            cor = abs(x_rec + A' * res);
            [a, b] = maxk(cor, idx);
            
            cor_s = x_rec + A' * res;
            x_rec(b, :) = cor_s(b, :);
            res = b_arr(:, idx) - A * x_rec;
            count = count + 1;
        end
        x_reconstructed(:, idx) = x_rec;
    end
end 