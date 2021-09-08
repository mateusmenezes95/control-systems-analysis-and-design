clear all, close all, clc
addpath(genpath("../../lib")) % Add lib path to Octave script file search paths

load_common_functions

% ===============================================================================
% Simulation parameters
% ===============================================================================

sim.time_step = 0.01;
sim.time = 50;
sim.sim_time = get_sim_time (sim.time_step, sim.time);

% ===============================================================================
% Reference parameters
% ===============================================================================

reference.start_time = 2;
reference.end_time = inf;
reference.amplitude = 1;

% ===============================================================================
% Disturbance parameters
% ===============================================================================

input_disturbance.start_time = 15;
input_disturbance.end_time = inf;
input_disturbance.amplitude = -0.2;

output_disturbance.start_time = 25;
output_disturbance.end_time = inf;
output_disturbance.amplitude = -0.2;

% ===============================================================================
% Transfer Functions Definions
% ===============================================================================

K = 0.5;
G = 2 / s;
C = K;

T = minreal((K * C * G) / (1 + (K * C * G)))
input_disturbance.tf = minreal(1 / (1 + (C * G)))
output_disturbance.tf = minreal(G / (1 + (C * G)))

reference.signal = get_step_signal(sim.sim_time, reference.start_time, reference.end_time, reference.amplitude);
input_disturbance.signal = get_step_signal(sim.sim_time, input_disturbance.start_time, input_disturbance.end_time, input_disturbance.amplitude);
output_disturbance.signal = get_step_signal(sim.sim_time, output_disturbance.start_time, output_disturbance.end_time, output_disturbance.amplitude);

[reference.y, reference.t, reference.x] = lsim(T, reference.signal, sim.time);
[input_disturbance.y, input_disturbance.t, input_disturbance.x] = lsim(input_disturbance.tf, input_disturbance.signal, sim.time);
[output_disturbance.y, output_disturbance.t, output_disturbance.x] = lsim(output_disturbance.tf, output_disturbance.signal, sim.time);

y = reference.y + input_disturbance.y + output_disturbance.y;
t = reference.t;

nfont=15;   % define tamanho da fonte de texto
nlinha=3;   % define espessura da linha

figure(1)

subplot(2,1,1)
plot(t,reference.signal,'--r','linewidth',nlinha-0.5)
hold on
plot(t,y,'b','linewidth',nlinha)
grid on
axis([0 t(end) min(y)-0.2  max(y)+0.2])
set(gca,'fontsize',nfont)

% Define as entradas de texto e ajusta o tamanho da fonte
hx=xlabel('Tempo (s)');
hy=ylabel('Saída');
%hl=legend('r(t)', 'y(t)');
%set(hl,'fontsize',nfont)
%set(hl,'location','southeast')
set(hx,'fontsize',nfont)
set(hy,'fontsize',nfont)

