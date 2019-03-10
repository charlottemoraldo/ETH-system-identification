function [p1_G_y,p1_omega,p1_best_gamma,p1_G_y_best_hann] = HS2018_SysID_midterm_p1_15819790(p1_u,p1_y,p1_z)
close all
clc
%% --------------------------------------------------------------------- %%
%% System Identification Midterm Problem 1
%% --------------------------------------------------------------------- %%
% Inputs:
%       > "p1_u" given input signal u(k)
%       > "p1_y" system output y(k) in response to p1_u for one realization
%       of noise v(k) = v1(k)
%       > "p1_z" system output z(k) in response to p1_u for one realization
%       of noise v(k) = v2(k)
% Outputs:
%       > "p1_G_y" the unsmoothed plant estimate calculated form the output
%       y(k) at the frequencies p1_omega
%       > "p1_omega" at least 30 equally-spaced frequencies in the range
%       [0.01,0.8] [rad/s] where the plant estimates p1_g_y and
%       p1_g_y_best_hann are provided
%       > "p1_best_gamma" the optimal choice of the Hann window width based
%       on a cross-validation between the two experiments
%       > "p1_G_y_best_hann" the plant estimate corresponding to the output
%       y(k), smoothed using the optimal choice of Hann window width
% ----------------------------------------------------------------------- %
%% My details
disp(' ');
disp('    ###########################################################')
disp('    #         SYSTEM IDENTIFICATION MIDTERM PROBLEM 1         #')
disp('    #              Student Name: Charlotte Moraldo            #')
disp('    #                Legi Number: 15-819-790                  #')
disp('    ###########################################################')
disp(' ');


%% --------------------------------------------------------------------- %%
%% PART 1 
%% --------------------------------------------------------------------- %%
disp(' ');
disp('#####################################################################')
disp('#                             PART 1:                               #')
disp('# Estimating the transfer function G_y in order to minimize the MSE #')
disp('#####################################################################')
disp(' ');


%% 1.1 Initialization

disp('### 1.1 Initialization ###############################################')
disp(' ');

%Hard code data from the midterm instructions
Ts = 2;
window_min = 0.01;
window_max = 0.8;
disp('I begin by hardcoding the following values, given in the midterm instructions:')
disp('    > Ts = 2s, the sampling time')
disp('    > [window_min, window_max] = [0.01, 0.8], the desired frequency range')
disp(' ');

%Find N the number of samples per period and N_per the number of periods
N = seqperiod(p1_u);
N_per_notround = size(p1_u,1)/seqperiod(p1_u);
N_per = floor(size(p1_u,1)/seqperiod(p1_u));
pts_discarded = size(p1_u,1) - N*N_per;
explanation3 = [...
    'After noticing that the signals u, y, and z are periodic, I compute their \n' ...
    'periodicity characteristics: \n' ...
    '    > The number of samples per period is obtained with the matlab \n'...
    '      command seqperiod(p1_u), where p1_u is given. I obtain: N = ',num2str(N),'.\n' ...
    '    > The number of periods is obtained through the division of the size \n'...
    '      of p1_u with the number of samples per period N: N_per = ',num2str(N_per_notround),'\n'...
    '      However, as we want a finite number of periods, we round this down to: \n' ...
    '      N_per = ',num2str(N_per),'\n'...
    '    > The number of points neglected is therefore equal to: \n'...
    '      1200 - N*N_per = ',num2str(pts_discarded),', where 1200 is the length of the input signals\n'...
    ];
fprintf(explanation3);
disp(' ');

%Extract the desired frequencies
omega = (2*pi/N)*[0:N-1]';
idx = find(omega>window_min*Ts & omega<window_max*Ts);
M = size(idx,1);
explanation4 = [...
    'To compute the vector p1_omega, I firstly compute a vector containing \n'...
    'N frequencies equally spread from 0 to 2*pi. \n' ...
    'I extract from this vector the frequencies which are contained in the \n'...
    'range Ts*[window_min, window_max] = Ts*[0.01, 0.8] = [0.02 1.6] rad/sample.\n'...
    'It is important not to forget the multiplication of the window by the\n'...
    'sampling time in order for the results to be coherent.\n'...
    ];
fprintf(explanation4);
disp(' ');
disp(' ');


%% 1.2 Estimating G_y

%1.2.1 Discarding the first 180 points, then using 4 periods
%Period 1
U1 = fft(p1_u(pts_discarded+1:N+pts_discarded)); 
Y1 = fft(p1_y(pts_discarded+1:N+pts_discarded));
G_etfe_y1 = Y1./U1;
%Period 2
U2 = fft(p1_u(N+pts_discarded+1:2*N+pts_discarded));
Y2 = fft(p1_y(N+pts_discarded+1:2*N+pts_discarded));
G_etfe_y2 = Y2./U2;
%Period 3
U3 = fft(p1_u(2*N+pts_discarded+1:3*N+pts_discarded));
Y3 = fft(p1_y(2*N+pts_discarded+1:3*N+pts_discarded));
G_etfe_y3 = Y3./U3;
%Period 4
U4 = fft(p1_u(3*N+pts_discarded+1:4*N+pts_discarded));
Y4 = fft(p1_y(3*N+pts_discarded+1:4*N+pts_discarded));
G_etfe_y4 = Y4./U4;
%Average over 4 periods
G_y_est1 = (G_etfe_y1 + G_etfe_y2 + G_etfe_y3 + G_etfe_y4)/4;
length1 = N*4;

%1.2.2 Discarding the first 180 points and the first period
G_y_est2 = (G_etfe_y2 + G_etfe_y3 + G_etfe_y4)/3;

disp('### 1.2 Estimating G_y ##############################################')
disp(' ');

explanation5 = ['In order to minimize the MSE, I must try to minimize the\n'...
               'effect of the noise and of the transient. When measuring a signal,\n'...
               'there is a high transient at the beginning due to inputs that might have\n'...
               'happened before the measurements and that havent been taken\n'...
               'into account. This transient spreads throughout the signal,\n'...
               'but decreases over the samples. \n'...
               '\n'...
               'As we have the advantage of having a periodic signal, I decided\n'...
               'to split the data into 4 records of one period each (N=255 samples),\n'...
               'which we can then average to get rid of the noise as much as possible.\n'...
               'There are several ways of doing so:\n'...
               ' \n'...
               '  1. At first thought, I would simply neglect the first period (samples 1 - ',num2str(N),')\n'...
               '     to reduce the transient, then average over the next 3 periods. \n'...
               '     This method however has the disadvantage that we neglect the last ',num2str(pts_discarded),'\n'...
               '     samples, which are useful data.\n'...
               '\n'...
               '  2. Instead, I can use the fact that the first samples are affected by \n'...
               '     the transient in order to neglect the first ',num2str(pts_discarded),' points.\n'...
               '     I would then still have 4 entire periods left of useful data. \n'...
               '     Using this, I still have the choice between: \n'...
               '        2.1 Neglecting the first ',num2str(pts_discarded),' samples as well as the first\n'...
               '            period in order to make sure that the effect of the transient\n'...
               '            is minimized. The average would then be computed over 3 periods. \n'...
               '        2.2 Only neglecting the first ',num2str(pts_discarded),' samples and average \n'...
               '            over 4 periods. \n'...
               '\n'...
               'A trade-off must be found between methods 2.1 and 2.2. While 2.1 could \n'...
               'be better in case the transient is very long, it also has the\n'...
               'the disadvantage of neglecting a lot of data points. In order to \n'...
               'keep as much useful data as possible, I therefore decided to keep \n'...
               'method 2.2, where we can make use of 4 entire periods of data, \n'...
               'corresponding to N*N_per = ',num2str(length1),' samples.\n'...
               'The averaged estimation obtained is returned in the variable p1_G_y.\n'...
    ];
fprintf(explanation5);
disp(' ');

%% 1.4 Plots
figure(1);
loglog(1/Ts*omega(idx), abs(G_etfe_y4(idx)),'linewidth',1); hold on; grid on;
loglog(1/Ts*omega(idx), abs(G_y_est1(idx)),'linewidth',1); hold on;
legend('Unaveraged ETFE of G_y', 'Averaged ETFE of G_y with 4 periods','location','southwest')
xlabel('log(w) [rad/s]');
ylabel('log|(G_y(w)|');
title('Estimate of transfer function G_y');
ax = gca;
ax.FontSize = 12;
axis([1/Ts*omega(idx(1)) 1/Ts*omega(idx(end)) 0 10])

explanation6 = ['In order to visualize the results of my implementation, the \n'...
    'plot (figure 1) shows: \n'... 
    '    > the unaveraged ETFE of G_y (taken over only data record)\n' ...
    '    > the averaged ETFE of G_y (using method 2.2, so over 4 periods) \n' ...
    ];
fprintf(explanation6);
disp(' ');

%% 1.5 Return outputs
p1_omega = omega(idx);
p1_G_y = G_y_est1;


%% --------------------------------------------------------------------- %%
%% PART 2 
%% --------------------------------------------------------------------- %%
disp(' ');
disp('#####################################################################')
disp('#                             PART 2:                               #')
disp('# Smoothing the estimate p1_G_y with a frequency-domain Hann window #')
disp('#####################################################################')
disp(' ');

%% 2.1 Estimage G_z similarly as in 1.3 and 1.4

% 2.1.1 Discarding the first 180 points and averaging over 4 periods
%Period 1
U1 = fft(p1_u(pts_discarded+1:N+pts_discarded)); 
Z1 = fft(p1_z(pts_discarded+1:N+pts_discarded));
G_etfe_z1 = Z1./U1;
%Period 2
U2 = fft(p1_u(N+pts_discarded+1:2*N+pts_discarded));
Z2 = fft(p1_z(N+pts_discarded+1:2*N+pts_discarded));
G_etfe_z2 = Z2./U2;
%Period 3
U3 = fft(p1_u(2*N+pts_discarded+1:3*N+pts_discarded));
Z3 = fft(p1_z(2*N+pts_discarded+1:3*N+pts_discarded));
G_etfe_z3 = Z3./U3;
%Period 4
U4 = fft(p1_u(3*N+pts_discarded+1:4*N+pts_discarded));
Z4 = fft(p1_z(3*N+pts_discarded+1:4*N+pts_discarded));
G_etfe_z4 = Z4./U4;
%Average over 4 periods
G_z_est1 = (G_etfe_z1 + G_etfe_z2 + G_etfe_z3 + G_etfe_z4)/4 ;

%2.2.2 Discarding the first 180 points and the first period
G_z_est2 = (G_etfe_z2 + G_etfe_z3 + G_etfe_z4)/3;

figure(2);
loglog(1/Ts*omega(idx), abs(G_etfe_z4(idx)),'linewidth',1); hold on; grid on;
loglog(1/Ts*omega(idx), abs(G_z_est1(idx)),'linewidth',1); hold on;
legend('Unaveraged ETFE of G_z', 'Averaged ETFE of G_z with 4 periods','location','southwest')
xlabel('log(w) [rad/s]');
ylabel('log|(G_z(w)|');
title('Estimate of transfer function G_z');
ax = gca;
ax.FontSize = 12;
axis([1/Ts*omega(idx(1)) 1/Ts*omega(idx(end)) 0 10])

disp('### 2.1 Estimating G_z ##############################################')
disp(' ');
explanation8 = ...
    ['In order to compute the estimate of G_z, the same method as part 1 is applied.\n'...
    'The first ',num2str(pts_discarded),' samples are discarded, then the ETFE is computed over \n'...
    '4 records each of size N = 255 (1 period), and they are then averaged.\n'...
    'Similarly as in part 1, figure 2 shows: \n'...
    '    > the unaveraged ETFE of G_z (taken over only data record)\n' ...
    '    > the averaged ETFE of G_z (over 4 periods) \n' ...
    ];
fprintf(explanation8);
disp(' ');

%% 2.2 Frequency-domain Hann windowing
figure(3); grid on;
plot_index = 1;
mse_tab = [];
window_sizes = [10,20,50,100,200];
best_Gs = 0*G_y_est1;
best_mse = 0;
for gamma = window_sizes
    [om, Wg_freq] = WfHann(gamma,N);
    Gs = 0*G_y_est1;
    zidx = find(om==0);
    om = [om(zidx:N);om(1:zidx-1)];
    Wg_freq = [Wg_freq(zidx:N) Wg_freq(1:zidx-1)];
    a = U1.*conj(U1);
    
    for wn=1:N
        norm=0;
        for xi=1:N
            widx = mod(xi-wn,N)+1;
            Gs(wn) = Gs(wn) + Wg_freq(widx)*G_y_est1(xi)*a(xi);
            norm = norm + Wg_freq(widx)*a(xi);
        end
        Gs(wn) = Gs(wn)/norm;
    end
    
    %Mean square error
    new_mse = 0; 
    Err = abs(Gs(idx)-G_z_est1(idx));
    for wn=1:size(idx,1)
        new_mse = new_mse + Err(wn)^2;
    end
    new_mse = new_mse/size(idx,1);
    mse_tab = [mse_tab new_mse];

    if plot_index == 1
        best_Gs = Gs;
        best_mse = new_mse;
    end
    if new_mse <= best_mse
        best_Gs = Gs;
        best_mse = new_mse;
    end
    
    %Plot smoothed estimate of G_y and unsmoothed G_z
    loglog(1/Ts*omega(idx),abs(Gs(idx)),'linewidth',1); hold on;  
    plot_index = plot_index+1;
end
[min_mse, min_index] = min(mse_tab);

loglog(1/Ts*omega(idx), abs(G_z_est1(idx)),'linewidth',1);
legend(['Smoothed estimate of G_y for gamma = ',num2str(window_sizes(1))],...
       ['Smoothed estimate of G_y for gamma = ',num2str(window_sizes(2))],...
       ['Smoothed estimate of G_y for gamma = ',num2str(window_sizes(3))],...
       ['Smoothed estimate of G_y for gamma = ',num2str(window_sizes(4))],...
       ['Smoothed estimate of G_y for gamma = ',num2str(window_sizes(5))],...
        'Unsmoothed estimate of G_z',...
        'Location','southwest');  
title('Unsmoothed estimate of G_z and smoothed estimate of G_y')
xlabel('log(w) [rad/s]');
ylabel('log|(G(w)|');
ax = gca;
ax.FontSize = 12;
axis([1/Ts*omega(idx(1)) 1/Ts*omega(idx(end)) 0 10])

disp(' ');
disp('### 2.2 Frequency-domain Hann window ################################')
disp(' ');
explanation9 = ...
    ['Using the given function WfHann, I now compute the Hann window, for \n' ...
    'gamma = [10,20,50,100,200]. I then apply it on the averaged estimate of G_y, in\n'...
    'order to smooth it, as discussed in class (slides 4.16 and 4.22) Figure 3  \n'...
    'shows the obtained windowed estimate for each gamma, as well as the unsmoothed\n'...
    'estimate of G_z. For each gamma, I also computed the mean-square error between \n'...
    'the y data estimate and the estimate based on the z data. The results obtained are:\n'...
    '    > gamma = ',num2str(window_sizes(1)),': MSE = ',num2str(mse_tab(1)),'\n'...
    '    > gamma = ',num2str(window_sizes(2)),': MSE = ',num2str(mse_tab(2)),'\n'...
    '    > gamma = ',num2str(window_sizes(3)),': MSE = ',num2str(mse_tab(3)),'\n'...
    '    > gamma = ',num2str(window_sizes(4)),': MSE = ',num2str(mse_tab(4)),'\n'...
    '    > gamma = ',num2str(window_sizes(5)),': MSE = ',num2str(mse_tab(5)),'\n'...
    'We can see that the minimum MSE is MSE = ',num2str(min_mse),', obtained with\n' ...
    'gamma = ',num2str(window_sizes(min_index)),', which is returned in the variable p1_best_gamma.\n' ...
    '\n' ...
    'Figure 4 illustrates the evolution of the MSE with respect to gamma. The\n' ...
    'optimal gamma for the minimum MSE is indicated in red.\n' ...
    ];
fprintf(explanation9);
disp(' ');

%% 2.3 Plot best estimates
figure(4);
plot(window_sizes,mse_tab,'b+-','linewidth',1); hold on; grid on;
plot(window_sizes(min_index),min_mse,'r*')
title('Evolution of the MSE with the Hann window width')
xlabel('Window size');
ylabel('Mean Square Error');
ax = gca;
ax.FontSize = 12;
axis([0 window_sizes(end) 0 mse_tab(1)])

%% 2.4 Return outputs
p1_best_gamma = window_sizes(min_index);
p1_G_y_best_hann = best_Gs;

%% End of function
end  