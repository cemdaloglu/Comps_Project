clear variables;
close all;
cvx_quiet TRUE;
%% 
n = 2^6;
m = 2^5;
A_random = randn([m, n]);
A_random = normc(A_random);

sparse_vector = zeros([n, 1]);

s = 6;
nnz_values = randn([s, 1]);
index = randsample(1:length(sparse_vector), s, true);
sparse_vector(index, 1) = nnz_values;

xx = A_random * sparse_vector;
b_random = awgn(xx, 10);
%% 

p = 0.1;
q = 0.1;
r = 1;
gamma = 0.1;

x0 = zeros(n, 1);
x = x0;
y = x0;

% gradient
GradF = @(x) A_random' * (A_random * x - b_random);
R = @(x) x0 * x; % !! FIX ME  !!
FBO = @(y) prox(gamma, R(x), y - gamma * GradF(y));

t = 1;

its = 1;
maxits = 10;
while(its<maxits)
    
    x_old = x;
    
    x = FBO(y);
    
    t_old = t;
    t = (p + sqrt(q + r * t_old^2)) / 2;
    a = (t_old - 1) / t;
    
    y = x + a * (x - x_old);
    
    res = norm(x_old-x, 'fro');
    
    
    if res / prod(n) < tol
        break;
    end
    
    its = its + 1;
    
end
%%
function [x] = prox(gamma, R, input)
    cvx begin
        variable x_hat(n, 1)
        minimize(gamma * R + 1/2 * norm(x_hat - input, 2)^2)
    cvx_end
    x = x_hat;
end