run common_parameters_script.m

observable_poles = filt(zpk(desirable_z_poles(3:4), [0 0], 1,sampling_period));
desirable_closed_loop_poles = controllable_poles*observable_poles;
desirable_coeficcients = get(desirable_closed_loop_poles, 'num'){1};

[a1, a2] = deal(a{1}(2), a{1}(3));
[b1, b2] = deal(b{1}(2), b{1}(3));

A = [
      b1  0   0   -1;
      b2  b1  0   1-a1;
      0   b2  b1  a1-a2;
      0   0   b2   a2
    ]

B = [
      desirable_coeficcients(2)-a1+1;
      desirable_coeficcients(3)+a1-a2;
      desirable_coeficcients(4)+a2;
      desirable_coeficcients(5)
    ]

x = inv(A)*B
[s0, s1, s2, r] = deal(x(1), x(2), x(3), x(4));

cz = filt([s0 s1 s2], [1 -r-1 r], sampling_period);

sz = filt([s0 s1 s2], 1, sampling_period)
rz = filt([1 -r-1 r], 1, sampling_period)
to = 1/dcgain(bz/controllable_poles)
tz = to*observable_poles
fz = minreal(tz/sz, minreal_precision)

pole_allocation_controller = cz;
reference_filter = fz;
save ./controllers/pole_allocation_controller.mat pole_allocation_controller reference_filter;

closed_loop = minreal((tz*bz/(rz*az+sz*bz)), minreal_precision)

[U, Y, E, R] = simulate_discrete_sys(sim.time, dt, integration_step_ratio,
                                  gn_of_s, cz, fz,
                                  reference.signal,
                                  output_disturbance.signal,
                                  input_disturbance.signal);

plot_response_and_control_signals (sim.time, figure_num = 1,
                                   reference.signal,
                                   Y, '$y(t)$',
                                   E, '$e(t)$',
                                   U, '$u(t)$',
                                   'b')
