clear variables;
close all;
cvx_quiet True;
%% 
S = load('data3_440.mat');
A_random = S.trainimages;
b = S.testimages;
b_random = b(:, 1);
%% 

p = 0.1;
q = 0.5;
r = 4;
%gamma = 0.1;
gamma = 1 /norm(A_random)^2;

x0 = zeros(size(A_random, 2), 1);
x = x0;
y = x0;

% 4.1 variables
e = randn([size(A_random, 1), 1]);
e = normc(e);
w = [x0; e];

idtt = speye(size(A_random, 1));
B = [A_random, idtt];

% w_hat = modified_fista(B, b(:, 1), w, gamma, p, q, r);
w_hat = fistamod(B, b(:, 1), w, gamma, p, q, r);
x_hat = w_hat(1:440, 1);
x_hat(x_hat < 0.001) = 0;
%%
function x_sparse = modified_fista(A_random, b_random, x, gamma, p, q, r)
    % gradient
    GradF = @(x) A_random' * (A_random * x - b_random);

    % Proxy
    % lambda = 10;
    % R = @(x) lambda .* norm(x, 1);
    % FBO = @(x, y) prox(y - gamma * GradF(y), gamma, R(x));


    FBO = @(y) max(y - gamma * GradF(y), gamma);
    
    tol = 10^(-6);
    t = 1;
    y = x;
    
    its = 1;
    maxits = 10;
    while(its<maxits)
        
        x_old = x;
        
        x = FBO(y);
        
        t_old = t;
        t = (p + sqrt(q + r * t_old^2)) / 2;
        a = (t_old - 1) / t;
        
        y = x + a * (x - x_old);
        
        res = norm(x_old - x, 'fro');
        
        
        if res / prod(440) < tol
            break;
        end
        its = its + 1; 
    end
    x_sparse = x;
end
%%
function x_sparse = fistamod(A_random, b_random, x_sprs, gamma, p, q, r)
    % gradient
    GradF = @(x) A_random' * (A_random * x - b_random);

    % Proxy
    lambda = 10;
    R = @(x) lambda .* norm(x, 1);
    FBO = @(x, y) prox(y - gamma * GradF(y), gamma, R(x));
    
    tol = 10^(-6);
    t = 1;
    x = x_sprs;
    y = x;
    
    its = 1;
    maxits = 100;
    while(its<maxits)
        
        x_old = x;
        
        x = FBO(x, y);
        
        t_old = t;
        t = (p + sqrt(q + r * t_old^2)) / 2;
        a = (t_old - 1) / t;
        
        y = x + a * (x - x_old);
        
        res = norm(x_old - x, 'fro');
        
        
        if res / 440 < tol
            break;
        end
        its = its + 1; 
    end
    x_sparse = x;
end
%%
function x = prox(input, gamma, R)
    cvx_begin
        variable x_hat1(40440, 1)
        minimize(gamma * R + 1/2 * pow_pos(norm(x_hat1 - input, 2), 2))
    cvx_end
    x = x_hat1;
end