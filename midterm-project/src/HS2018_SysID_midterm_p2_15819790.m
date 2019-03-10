function [p2_T, p2_u] = HS2018_SysID_midterm_p2_15819790(p2_G_mag, p2_omega, p2_y_step, p2_t_step)
%% --------------------------------------------------------------------- %%
%% System Identification Midterm Problem 2
%% --------------------------------------------------------------------- %%
% Inputs:
%       > "p2_G_mag" contains the magnitudes of G_approx corresponding to
%       the frequencies p2_omega in [rad/s]
%       > "p2_y_step" contains the step response values corresponding to
%       the times in p2_t_step
% Outputs:
%       > "p2_T" sampling time ([sec]) chosen for the experiment, as scalar
%       value
%       > "p2_u" chosen input signal as Nx1 vector
% ----------------------------------------------------------------------- %
%% My details
disp(' ');
disp('    ###########################################################')
disp('    #         SYSTEM IDENTIFICATION MIDTERM PROBLEM 2         #')
disp('    #              Student Name: Charlotte Moraldo            #')
disp('    #                Legi Number: 15-819-790                  #')
disp('    ###########################################################')
disp(' ');

%% Question A
disp(' ');
disp('#####################################################################')
disp('#         Question A: Characteristics of the input signal           #')
disp('#####################################################################')
disp(' ');
disp('The input signal p2_u is a PRBS signal of 21 periods, where each period')
disp('contains 1023 data points sampled every 74ms.')
disp(' ');

disp(' ');
disp('#####################################################################')
disp('#       Question B: Estimating G from the experimental data         #')
disp('#####################################################################')
disp(' ');
disp('In order to estimate G from the experimental data, I would do similarly')
disp('as what was done in the exercise 1:')
disp('  > Take advantage of the signal periodicity and split the input and')
disp('    output signals in 21 records of 1023 samples each (which')
disp('    corresponds to the length of one period).')
disp('  > Take the fft of the input and output signal over each record,')
disp('    then compute the ETFE and average them to obtain the estimate.')
disp('  > In order to further reduce the MSE, finish by smoothing the estimate')
disp('    by applying a frequency-domain Hann window. The window width gamma')
disp('    has to be optimized through testing in order to find the minimal MSE.')
disp(' ');

disp(' ');
disp('#####################################################################')
disp('#          Question C: Satisfying the given requirements            #')
disp('#####################################################################')
disp(' ');
%% Requirement 1
% The Nyquist frequency, and the highest frequency in your estimate, is 
% greater than the frequency at which the plant has rolled off to less than 
% 5% of its DC gain margin
dc_gain = 0.05 * p2_G_mag(1);
idx_dc = find(p2_G_mag < dc_gain);
idx_dc = idx_dc(1,1);
omega_dc = p2_omega(idx_dc);

% Justification:
disp('--------------------------- REQUIREMENT 1 --------------------------- ') 
disp('After using p2_G_mag to find that the DC gain is equal to 1.0894, I')
disp('use the command find() to get the highest frequency in p2_omega for')
disp('which the gain of p2_G_mag is lower than 5 percent of this DC gain,')
disp('which gives me a frequency of 41.969 rad/s.')
disp('The Nyquist frequency pi/p2_T therefore has to be greater than 41.969 rad/s,')
disp('from which we can deduce that the sampling time p2_T has to be smaller')
disp('than 0.0749 seconds, or 74.9 ms.')
disp(' ');
disp(' ')

%% Requirement 2
% p2_T is an integer number of millisecond
T_max = floor(pi/omega_dc*1000)/1000;

% Justification:
disp('--------------------------- REQUIREMENT 2 --------------------------- ') 
disp('From requirement 1, I derived that the sampling time has to be smaller')
disp('than 74.9 ms. Requirement 2 states that an integer number of milliseconds')
disp('is desired; furthermore, it is in our best interest to keep the largest')
disp('sampling time possible in order to optimize the data storage, therefore')
disp('I round down the result from requirement 1 and get a sampling time of')
disp('p2_T = 74ms.')
disp(' ');
disp(' ')

%% Requirement 3
% G_hat is an unbiased estimate of G. Justification:
disp('--------------------------- REQUIREMENT 3 --------------------------- ') 
disp('From what was seen in class (slide 3.15), I know that for periodic inputs ');
disp('(where there is an integer number of periods), the ETFE is unbiased.')
disp('Therefore, I decided to design my input as a PRBS signal, as it has');
disp('interesting identification properties, for example a constant PSD for');
disp('finite length signals.');
disp(' ')
disp(' ')

%% Requirement 4
% Frequency resolution
idx_peak = find(p2_G_mag == max(p2_G_mag));
omega_peak = p2_omega(idx_peak);
N_max = round(4*2*pi/(T_max*omega_peak));
x = ceil(log2(N_max + 1));
N = 2.^x -1;

%Justification:
disp('--------------------------- REQUIREMENT 4 --------------------------- ') 
disp('After using the command find() to get the low-frequency peak of p2_G_mag')
disp('(equal to 0.5549 rad/s), I can use the definition of the frequency')
disp('resolution (slide 1.32) to write: 4 * 2*pi/(N*T) < 0.5549 (where T=74ms, N')
disp('is the number of samples per period, and the factor 4 represents the 4 non-zero')
disp('frequencies we want to have below 0.5549) -- this inequality yields: N > 613.')
disp(' ')
disp('Knowing that the period length of a PRBS must be at most N = 2^x - 1,')
disp('we find that for N > 613 we have x > 9.26, which we round up to 10 since')
disp('x should be an integer and N has to be strictly bigger than 613 -- this')
disp('finally yields that every period of the signal will contain:')
disp('N = 2^10-1 =  1023 samples.')
disp(' ')
disp(' ')

%% Requirement 5
% Variance
var_input = 0.025;
var_noise = 0.5;
N_per = var_noise/var_input+1;

%Justification:
disp('--------------------------- REQUIREMENT 5 --------------------------- ') 
disp('I know from the class (slide 4.4) that when splitting the data into');
disp('N_per segments of length N = 1023 (one period) and then averaging their');
disp('ETFE, the variance decays proportionally to 1/N_per.')
disp(' ')
disp('Furthermore, the signal y = Gu + v has the same noise distribution as')
disp('v, which means that the variance of y is equal to the variance of v,')
disp('thus:')
disp('   var_noise / N_per < var_error --> N_per > 0.5/0.025 = 20 periods')
disp('and we therefore choose a number of period equal to 21 in order to')
disp('satisfy this inequality (as we want the variance to be strictly below')
disp('0.025).')

%% Build input signal
%PRBS signal of N_per = 21 periods, each of N = 1023 samples
u = idinput([N 1 1],'prbs',[0 1],[-2.5 2.5]);
u_per = [];
for i=1:N_per
    u_per = [u_per; u];
end

%% Return outputs
p2_T = T_max
p2_u = u_per;
size(p2_u)

%% End of function
end 

