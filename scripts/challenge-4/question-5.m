run common_parameters_script.m

load ./controllers/pole_allocation_controller.mat
load ./controllers/imc_controller.mat

[U, Y, E, R] = simulate_discrete_sys(sim.time, dt, integration_step_ratio,
                                  gn_of_s,
                                  pole_allocation_controller,
                                  reference_filter,
                                  reference.signal,
                                  output_disturbance.signal,
                                  input_disturbance.signal);

plot_response_and_control_signals (sim.time, figure_num = 1,
                                   reference.signal,
                                   Y, 'Alocação de Polos',
                                   E, 'Alocação de Polos',
                                   U, 'Alocação de Polos',
                                   'b')

hold on

[U, Y, E, R] = simulate_discrete_sys(sim.time, dt, integration_step_ratio,
                                  gn_of_s,
                                  imc_controller,
                                  filt(1),
                                  reference.signal,
                                  output_disturbance.signal,
                                  input_disturbance.signal);

plot_response_and_control_signals (sim.time, figure_num = 1,
                                   reference.signal,
                                   Y, 'IMC',
                                   E, 'IMC',
                                   U, 'IMC',
                                   'k')

subplot(3,1,1)
line([4 4], ylim(), 'color', 'k', 'linestyle', '--')
set(legend, 'string', get(legend,'string')(1:3))
set(legend,'orientation', 'horizontal')
set(legend,'location', 'southeast')

subplot(3,1,2)
set(legend,'orientation', 'horizontal')
set(legend,'location', 'northeast')

subplot(3,1,3)
set(legend,'orientation', 'horizontal')
set(legend,'location', 'southwest')
