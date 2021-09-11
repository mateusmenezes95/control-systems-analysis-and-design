clear all, close all, clc
addpath(genpath("../../lib")) % Add lib path to Octave script file search paths

run common_functions_script

% ===============================================================================
% Simulation parameters
% ===============================================================================

sim = get_sim_time (dt = 0.01, end_time = 50);

% ===============================================================================
% Input signals
% ===============================================================================

unit_step = get_signal (sim.time, start_time = 2, end_time = inf, amplitude = 1);
output_disturbance = get_signal (sim.time, start_time = 15, end_time = inf, amplitude = -0.2);
input_disturbance = get_signal (sim.time, start_time = 25, end_time = inf, amplitude = -0.2);

% ===============================================================================
% Transfer Functions Definions
% ===============================================================================

F = 1;
K = 0.5;
G = 2 / s;
C = K;

% ===============================================================================
% Main of the script
% ===============================================================================

y_to_r = get_output_to_reference_tf (G, C, F)
y_to_qy = get_output_to_output_disturbance_tf (G, C, F)
y_to_qu = get_output_to_input_disturbance_tf (G, C, F)

unit_step.response = lsim(y_to_r, unit_step.signal, sim.time);
output_disturbance.response = lsim(y_to_qy, output_disturbance.signal, sim.time);
input_disturbance.response = lsim(y_to_qu, input_disturbance.signal, sim.time);

control_loop_response = unit_step.response + output_disturbance.response + input_disturbance.response;

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

