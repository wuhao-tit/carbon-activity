clear;
clc;
%% activity of C calculate
W = 8389.31; %J/mol
R = 8.314; %J/(mol*K)
T = 1000; %K
% x = 0.05:0.001:0.012;
x = 0.05;
% syms x
AJ= 1 - exp(-W/(R*T));
delta = sqrt(1 - 2*(1+2*AJ)*x + (1+8*AJ)*x^2);
activity_C = 5 * log((1-2*x)/x) + 6*W/(R*T) + (38575-13.48*T)/(R*T)...
    + 6 * log((delta-1+3*x)/(delta+1-3*x));
%% DEF activity
DEF_delta = (0.5 / delta) * (-2 - 4*AJ + 2*x + 16.0*AJ*x);
DEF_activity_C = -((10/(1 - 2*x)) + (5/x)) + 6*((DEF_delta + 3)...
     /(delta-1+3*x) - (DEF_delta-3) / (delta+1-3*x));
%% diffusion coefficient calc
z = 12; % for austenite
Boltz = 1.380649e-23; %J/K
HH = 6.62607015e-34; %Js;
theta = x/(1-x);% c atoms / Fe atoms
activity_C = exp(activity_C);%activity_C transform from log to nature
DEF_activity_C = DEF_activity_C * activity_C;
DEF_activity_C = DEF_activity_C * 1 / (1+theta)^2;
niu = activity_C * (1 + z*(1+theta) / ((1 - (0.5*z+1)*theta +(0.25*z^2+0.5*z)*(1-AJ)*theta^2)))...
    + (1+theta) * DEF_activity_C;
Diff = Boltz*T/HH *exp(-21230/T)*exp(-31.84)*niu *1e-4; %C diffusion coefficient, m**/s
% %% calc avg_diffusion coefficient
% syms x_surf x_inter x_interval i xx
% %division the x from surf to inter, calc the Diffusion coefficient at each
% %x, then calc D_surf to D_inter
% D_avg =D_integral /(x_surf - x_inter);
% D_integral = (sum(D)+ (D(1)+D(number_of_xx)/2))*x_inter;
% %% calc wt%2mol_frac
% %% lamda_crit_spa, which calculated from supercooling rate.
% %   the calc from