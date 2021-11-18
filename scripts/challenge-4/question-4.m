run common_parameters_script.m

kf = 1/dcgain(bz/controllable_poles)
fr_of_z = kf*(bz/controllable_poles)

[fr_num_coefficients, fr_den_coefficients, ts] = filtdata(fr_of_z);
nfr_of_z = filt(fr_num_coefficients{1}, 1, ts)
dfr_of_z = filt(fr_den_coefficients{1}, 1, ts)

cz = (kf*az)/(dfr_of_z - nfr_of_z)

imc_controller = cz;
save ./controllers/imc_controller.mat imc_controller;

closed_loop = minreal(feedback(cz*gn_of_z, 1), minreal_precision)
y_to_qy = minreal(1/(1+(cz*gn_of_z)), minreal_precision)
y_to_qu = minreal(gn_of_z/(1+(cz*gn_of_z)), minreal_precision)

[U, Y, E, R] = simulate_discrete_sys(sim.time, dt, integration_step_ratio,
                                  gn_of_s, cz, filt(1),
                                  reference.signal,
                                  output_disturbance.signal,
                                  input_disturbance.signal);

plot_response_and_control_signals (sim.time, figure_num = 1,
                                   reference.signal,
                                   Y, '$y(t)$',
                                   E, '$e(t)$',
                                   U, '$u(t)$',
                                   'b')
