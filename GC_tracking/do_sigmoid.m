function [sig_fit, sig_gof] = do_sigmoid(time, speed)

%this function calculates a sigmoidal fit to twitching motility data

v_h = 1.5; %start point for high velocity mode
v_l = 0.5; %start point for loe velocity mode
T_gs = 6; %Time point of Global Switching in [min]
tau_gs = 0.5; %Global Switching time window in [min]

my_sigmoid = fittype('v_h - (v_h - v_l)./(1+exp((T_gs - x)./(0.17*tau_gs)))',...
   'coeff',{'v_h','v_l', 'T_gs', 'tau_gs'});

opts = fitoptions(my_sigmoid);
set(opts,'TolFun',1E-12, 'TolX', 1E-12, 'StartPoint', ...
    [v_h, v_l, T_gs, tau_gs], 'Lower', [0, 0, 0, 0], ...
    'Upper', [1.6, 2, 100, 200]);

[sig_fit, sig_gof] = fit(time, speed, my_sigmoid,opts);
