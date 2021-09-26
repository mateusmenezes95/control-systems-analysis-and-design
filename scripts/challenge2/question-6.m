addpath(genpath("./")) % Add lib path to Octave script file search paths
run common_parameters_scripts.m

% =============================================================================
% Transfer Functions Definions
% =============================================================================

disp('=======================================================================')
disp('Simulating plant response with input reference filter')
disp('=======================================================================')

model_num = 1
delay = l(model_num)
c = kc * ((s*ti + 1) / s*ti)
f = 1 / ((ti*s) + 1)
g = g = k(model_num) / (tau(model_num) * s + 1)

% =============================================================================
% Main of the script
% =============================================================================

% [u, y, err] = simulate_sys(sim.time, dt,
%                            g, delay, c, f,
%                            reference.signal,
%                            step_output_disturbance.signal,
%                            step_input_disturbance.signal);

% plot_response_and_control_signals (sim.time, figure_num = 1,
%                                   reference.signal,
%                                   y, '$y_{1}(t)$ com filtro de referência',
%                                   u, '$u_{1}(t)$ com filtro de referência', 'b-')

% disp('=======================================================================')
% disp('Simulating plant response with I+P controller')
% disp('=======================================================================')

f = 1
c = kc / (s*ti)

% [u, y, err] = simulate_sys(sim.time, dt,
%                            g, delay, c, f,
%                            reference.signal,
%                            step_output_disturbance.signal,
%                            step_input_disturbance.signal,
%                            controller_type = 'I+P');

% hold on

% plot_response_and_control_signals (sim.time, figure_num = 1, 'none',
%                                   y, '$y_{1}(t)$ com controlador I+P',
%                                   u, '$u_{1}(t)$ com controlador I+P', 'm-')

% disp('=======================================================================')
% disp('Simulating plant response with I+P controller and control outuput saturation')
% disp('=======================================================================')

% [u, y, err] = simulate_sys(sim.time, dt,
%                            g, delay, c, f,
%                            reference.signal,
%                            step_output_disturbance.signal,
%                            step_input_disturbance.signal,
%                            controller_type = 'I+P',
%                            saturation = [0 0.8]);

% hold on

% plot_response_and_control_signals (sim.time, figure_num = 1, 'none',
%                                   y, '$y_{1}(t)$ com saturação de $u(t)$',
%                                   u, '$u_{1}(t)$ com saturação', 'g-')

% disp('=======================================================================')
% disp('Simulating plant response with I+P controller and control outuput saturation ')
% disp('=======================================================================')

[u, y, err] = simulate_sys(sim.time, dt,
                           g, delay, c, f,
                           reference.signal,
                           step_output_disturbance.signal,
                           step_input_disturbance.signal,
                           controller_type = 'I+P',
                           saturation = [0 0.8],
                           anti_windup_tol = 10e-5);

hold on

plot_response_and_control_signals (sim.time, figure_num = 1, reference.signal,
                                  y, '$y_{1}(t)$ com anti windup',
                                  u, '$u_{1}(t)$ com anti windup', 'b-')

subplot(4,1,2)
plot_signal(sim.time, 'Tempo (s)',
            err, 'Erro $e(t)$',
            '-b', '$e(t)$')
