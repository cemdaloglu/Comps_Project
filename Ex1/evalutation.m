%% evaluation 

%%
titletoset = {'Averaged Normalized Reconstruction error as a ', 'function of sparsity s using RANDOM sensor matrix'};
legendtoinsert = {'l1', 'OMP', 'BT', 'MP', 'HTP', 'SP', 'CoSaMP', 'IHT'};
plot_rec_err_all({'l1_random', 'OMP_random', 'BT_random', 'MP_random', 'HTP_random', 'SP_random', 'CoSaMP_random','IHT_random'}, titletoset, legendtoinsert, [0,8])

%% 

titletoset = {'Averaged Normalized Reconstruction error as a', 'function of sparsity s using FOURIER sensor matrix'};
legendtoinsert = {'l1', 'OMP', 'MP', 'BT', 'HTP', 'SP', 'CoSaMP', 'IHT'};
plot_rec_err_all({'l1_fourier', 'OMP_fourier', 'MP_fourier','BT_fourier', 'HTP_fourier', 'SP_fourier', 'CoSaMP_fourier', 'IHT_fourier'}, titletoset, legendtoinsert, [0,1.5])


%% 
% 
function plot_rec_err_all(algorithm_cell_list, titletoset, legendtoinsert, ylims)
    figure; hold on

    x = load('x_true.mat');
    x_true = x.x_true; 


    for str = algorithm_cell_list
        % load x and reconstruction
        x_rec = load(strcat(str{1},'.mat')).x_recovered;

        x_x_rec_diff = abs(x_rec - x_true);
        rec_error = vecnorm(x_x_rec_diff, 2)./vecnorm(x_true, 2);
        rec_error_avg = mean(rec_error, 3);  
        
        s = 1:1:32;
        
        plot(s,rec_error_avg, '-x', 'LineWidth',1.5)
        title(titletoset)
        ylim(ylims)
        xlabel('sparsity s')
        ylabel('reconstruction error')
        
        
    end
    yline(10e-6) 
    legend(legendtoinsert, 'Location','northwest')
    hold off
    
end 







