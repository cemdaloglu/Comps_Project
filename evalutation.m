%% evaluation 

%% ﻿plot the averaged normalized recovery error 
% ||x−x_hat||_2=||x||_2 of the different reconstruction algorithms for each s.
% Declare a recovery successful if ||x − x_hat||2=||x||2 ≤ 10−6

x_x_rec_diff = x_recovered - x_true;
rec_error = vecnorm(x_x_rec_diff, 2)./vecnorm(x_true, 2)
rec_error_avg = mean(rec_error, 3)

%%
%size(x_recovered - x_true)
%size(vecnorm(x_recovered - x_true, 2)./vecnorm(x_true, 2))
size(mean(vecnorm(x_recovered - x_true, 2)./vecnorm(x_true, 2), 3))


%% Example for 100 repetitions of l1 with random matrix 
x = load('x_true.mat');
x_true = x.x_true; 

l1_x_rec = load(strcat('l1','.mat'));
l1_x_rec_random = l1_x_rec.x_recovered_random;

%plot_rec_err(x_true, l1_x_rec_random)

plot_rec_err_all({'l1', 'OMP'}, true, false)

%%
titletoset = {'Averaged Normalized Reconstruction error as a ', 'function of sparsity s using random sensor matrix'};
legendtoinsert = {'l1', 'OMP', 'BT', 'MP', 'HTP'};
plot_rec_err_all({'l1_random', 'OMP_random', 'BT_random', 'MP_random', 'HTP_random'}, titletoset, legendtoinsert)

%% 

titletoset = {'Averaged Normalized Reconstruction error as a', 'function of sparsity s using FOURIER sensor matrix'};
legendtoinsert = {'l1', 'MP', 'BT', 'MP', 'HTP'};
plot_rec_err_all({'l1_fourier', 'MP_fourier', 'BT_fourier', 'MP_random', 'HTP_random'}, titletoset, legendtoinsert)

%%
function plot_rec_err(x_true, x_rec_rep)
    x_x_rec_diff = x_rec_rep - x_true;
    rec_error = vecnorm(x_x_rec_diff, 2)./vecnorm(x_true, 2);
    rec_error_avg = mean(rec_error, 3);
    
    s = 1:1:32;
    figure
    %plot(s,mean_rec_error_random, s, mean_rec_error_random + 2 )
    plot(s,rec_error_avg, '-x' )
  

end 



%% 
% 
function plot_rec_err_all(algorithm_cell_list, titletoset, legendtoinsert)
    figure; hold on

    x = load('x_true.mat');
    x_true = x.x_true; 


    for str = algorithm_cell_list
        % load x and reconstruction
        x_rec = load(strcat(str{1},'.mat')).x_recovered;

        x_x_rec_diff = x_rec - x_true;
        rec_error = vecnorm(x_x_rec_diff, 2)./vecnorm(x_true, 2);
        rec_error_avg = mean(rec_error, 3);  

        s = 1:1:32;
        
        plot(s,rec_error_avg, '-x')
        title(titletoset)
        xlabel('sparsity s')
        ylabel('reconstruction error')
        
        
    end
    yline(10e-6) 
    legend(legendtoinsert)
    hold off
end 







