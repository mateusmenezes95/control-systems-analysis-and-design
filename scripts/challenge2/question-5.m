addpath(genpath("./")) % Add lib path to Octave script file search paths
run common_parameters_scripts.m

% =============================================================================
% Transfer Functions Definions
% =============================================================================

c = kc * ((s*ti + 1) / s*ti);
f = 1;

% =============================================================================
% Graph definitions
% =============================================================================

line_colors = {'-b;', '-g;', '-m;', '-k;'};

% =============================================================================
% Main of the script
% =============================================================================

subplot(2,1,1)
plot_signal(sim.time, 'Tempo (s)',
            reference.signal, 'Saída $y_{i}(t)$',
            '--r;$r(t)$;')
hold on

for i = 1:num_plants
    g = k(i) / (tau(i) * s + 1)
    [u, y, err] = simulate_sys(sim.time, dt,
                             g, l(i), c, f,
                            reference.signal,
                            step_output_disturbance.signal,
                            step_input_disturbance.signal);
    subplot(2,1,1)
    plot_signal(sim.time, 'Tempo (s)',
               y, 'Saída $y_{i}(t)$',
               [line_colors{i} '$y_{' num2str(i) '}(t)$;'])

    subplot(2,1,2)
    plot_signal(sim.time, 'Tempo (s)',
                u, 'Sinal $u_{i}(t)$',
                [line_colors{i} '$u_{' num2str(i) '}(t)$;'])
    
    hold on
endfor
