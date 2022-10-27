function [ x_reconstructed ] = SP_algorithm(A, b_arr)
    [m,n] = size(A);
    x_reconstructed = zeros(n, m);
    
    for idx = 1:m
        res = b_arr(:, idx);
        count = 0;
        while norm(res, 2) > 10^(-6) && count < m
            cor = abs(A' * res);
            [a, b] = maxk(cor, idx);
    
            cvx_begin
                variable x_rec1(idx, 1)
                minimize(norm(b_arr(:, idx) - A(:, b) * x_rec1, 2))
            cvx_end
              
            x_ = zeros(n, 1);
            x_(b, :) = x_rec1;
            [a2, b2] = maxk(abs(x_), idx);
    
    
            cvx_begin
                variable x_rec2(idx, 1)
                minimize(norm(b_arr(:, idx) - A(:, b) * x_rec2, 2))
            cvx_end
              
            x_sup = zeros(n, 1);
            x_sup(b, :) = x_rec2;
    
            res = b_arr(:, idx) - A * x_sup;
            count = count + 1;
        end
        x_reconstructed(:, idx) = x_sup;
    end

end 