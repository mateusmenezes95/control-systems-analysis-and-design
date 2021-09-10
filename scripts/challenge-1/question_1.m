clear all, close all, clc
addpath(genpath("../../lib")) % Add lib path to Octave script file search paths

load_common_functions

% ===============================================================================
% Simulation parameters
% ===============================================================================

sim = get_sim_time (dt = 0.01, end_time = 50);

% ===============================================================================
% Input signals
% ===============================================================================

unit_step = get_input_signal (sim.time, start_time = 2, end_time = inf, amplitude = 1);
input_disturbance = get_input_signal (sim.time, start_time = 15, end_time = inf, amplitude = -0.2);
output_disturbance = get_input_signal (sim.time, start_time = 25, end_time = inf, amplitude = -0.2);

% ===============================================================================
% Transfer Functions Definions
% ===============================================================================

K = 0.5;
G = 2 / s;
C = K;

Y_R = minreal((K * C * G) / (1 + (K * C * G)));  % Y(s)/R(s)
Y_R.inname = 'R(s)';
Y_R.outname = 'Y(s)'

Y_Qu = minreal(1 / (1 + (C * G)));               % Y(s)/Qu(s)
Y_Qu.inname = 'Qu(s)';
Y_Qu.outname = 'Y(s)'

Y_Qy = minreal(G / (1 + (C * G)));               % Y(s)/Qy(s)
Y_Qy.inname = 'Qy(s)';
Y_Qy.outname = 'Y(s)'

unit_step.response = lsim(Y_R, unit_step.signal, sim.time);
input_disturbance.response = lsim(Y_Qu, input_disturbance.signal, sim.time);
output_disturbance.response = lsim(Y_Qy, output_disturbance.signal, sim.time);

control_loop_response = unit_step.response + input_disturbance.response + output_disturbance.response;

nfont=15;   % define tamanho da fonte de texto
nlinha=2;   % define espessura da linha

figure(1)

subplot(2,1,1)
plot(sim.time, unit_step.signal, '--r', 'linewidth', nlinha-0.5)
hold on
plot(sim.time, control_loop_response, 'b', 'linewidth', nlinha)
grid on
axis([0 sim.time(end) min(control_loop_response)-0.2  max(control_loop_response)+0.2])
set(gca,'fontsize',nfont)

% Define as entradas de texto e ajusta o tamanho da fonte
hx=xlabel('Tempo (s)');
hy=ylabel('Saída');
%hl=legend('r(t)', 'y(t)');
%set(hl,'fontsize',nfont)
%set(hl,'location','southeast')
set(hx,'fontsize',nfont)
set(hy,'fontsize',nfont)

