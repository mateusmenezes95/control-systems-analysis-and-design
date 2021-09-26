printf("Loading commons functions...\n");

pkg load control
pkg load signal

s = tf("s");

function sim = get_sim_time (dt = 0.01, end_time = 100)
    sim.time_step = dt;
    sim.time = 0:dt:(end_time-dt);
    sim.time = transpose(sim.time);
    if (nargin == 0)
        printf("No arguments passed. Using default values instead\n")
    end
end

function [u, start_idx, end_idx]  = get_step_signal (sim_time, amplitude, start_time, end_time)
    sim_end_time = sim_time(length(sim_time));

    if end_time > sim_end_time
        disp(cstrcat(num2str(end_time), " > ", num2str(sim_end_time),
                     ". Then, end_time = ", num2str(sim_end_time)))
        end_time = sim_end_time;
    end

    u = zeros(length(sim_time), 1);
    start_idx = find(sim_time >= start_time)(1) + 1;
    end_idx = find(sim_time >= end_time)(1);
    u(start_idx:end_idx) = amplitude;
end

function ramp = get_ramp_signal (sim_time, slope, start_time, end_time)
    [ramp, start_idx, end_idx] = get_step_signal (sim_time, 1, start_time, end_time);
    ramp(start_idx:end_idx) = (slope * ramp(start_idx:end_idx)
                               .* (sim_time(start_idx:end_idx) - start_time));
end

function sine = get_sine_signal (sim_time, amplitude, natural_freq,
                                 start_time, end_time)
    [sine, start_idx, end_idx] = get_step_signal (sim_time, 1, start_time, end_time);
    sine(start_idx:end_idx) = amplitude * sin(natural_freq.*sim_time(start_idx:end_idx));
end

% function u = get_step_signal (sim_time, param = struct)
%     u = get_step_signal (sim_time, param.start_time,
%      param.end_time, param.amplitude)
% end

function input_signal = get_signal (sim_time, amplitude,
                                    start_time, end_time,
                                    signal = 'step')
    input_signal.start_time = start_time;
    input_signal.end_time = end_time;

    switch (signal)
        case 'step'
            input_signal.amplitude = amplitude;
            input_signal.signal = get_step_signal(sim_time, amplitude,
                                                  start_time, end_time);
        case 'ramp'
            input_signal.slope = amplitude;
            input_signal.signal = get_ramp_signal(sim_time, amplitude,
                                                  start_time, end_time);
        otherwise
            error('Choose invalid signal! Valids are step or ramp')
    endswitch
end

% Y(s)/R(s)
function y_to_r_tf = get_output_to_reference_tf (G, C, F)
    y_to_r_tf = minreal(F * ((C * G) / (1 + (C * G))));
    y_to_r_tf.inname = 'R(s)';
    y_to_r_tf.outname = 'Y(s)';
end

% Y(s)/Qu(s)
function y_to_qu_tf = get_output_to_input_disturbance_tf (G, C, F)
    y_to_qu_tf = minreal(G / (1 + (C * G)));
    y_to_qu_tf.inname = 'Qu(s)';
    y_to_qu_tf.outname = 'Y(s)';
end

% Y(s)/Qy(s)
function y_to_qy_tf = get_output_to_output_disturbance_tf (G, C, F)
    y_to_qy_tf = minreal(1 / (1 + (C * G)));
    y_to_qy_tf.inname = 'Qy(s)';
    y_to_qy_tf.outname = 'Y(s)';
end

% E(s)/R(s)
function e_to_r_tf = get_error_to_reference (G, C, F)
    e_to_r_tf = get_output_to_output_disturbance_tf(G, C, F);
    e_to_r_tf.inname = 'R(s)';
    e_to_r_tf.outname = 'E(s)';
end

 % U(s)/R(s)
function u_to_r_tf = get_controller_output_to_reference_tf (G, C, F)
    u_to_r_tf = minreal(F * (C / (1 + (C * G))));
    u_to_r_tf.inname = 'R(s)';
    u_to_r_tf.outname = 'U(s)';
end

 % U(s)/Qy(s)
function u_to_qu_tf = get_controller_output_to_output_disturbance_tf (G, C, F)
    u_to_qu_tf = minreal(((-1 * C) / (1 + (C * G))));
    u_to_qu_tf.inname = 'Qy(s)';
    u_to_qu_tf.outname = 'U(s)';
end

 % U(s)/Qu(s)
function u_to_qu_tf = get_controller_output_to_input_disturbance_tf (G, C, F)
    u_to_qu_tf = minreal(((-1 * C * G) / (1 + (C * G))));
    u_to_qu_tf.inname = 'Qu(s)';
    u_to_qu_tf.outname = 'U(s)';
end

function [xkplus1, y] = get_ss_output(x0, ss_matrix = ss(0,0,0,0), uk) 
    k = 1;
    xk = x0;
    xkplus1 = xk*(1+ss_matrix.a) + (ss_matrix.b*uk);
    y = ss_matrix.c*xkplus1 + (ss_matrix.d*uk);
end

% Discrete simulation based in SANTOS, T. L. M. code example 
function [U, Y, E] = simulate_sys(sim_time, dt, p, delay, c, f, r,
                                  qu = 0, qy = 0,
                                  controller_type = 'standard',
                                  saturation = [],
                                  anti_windup_tol = inf)
    if !isscalar(f)
        f = ss(f);
        f.a = dt*f.a;
        f.b = dt*f.b;
        xf = [zeros(size(f.a,1),1)];
    end

    p = ss(p);
    p.a = dt*p.a;
    p.b = dt*p.b;
    
    kc = get(c, 'num'){1};
    c = ss(c);
    c.a = dt*c.a;
    c.b = dt*c.b;

    ld = round(delay/dt);
    u_delay = zeros(1,ld+1);

    xc = [zeros(size(c.a,1),1)];
    xp = [zeros(size(p.a,1),1)];
    y = p.c * xp;

    u_set = 0;
    u_desirable = 0;
    x = zeros(1, length(sim_time));
    us = x;
    ud = x;
    aw = x;

    for k=1:length(sim_time)
        if isscalar(f)
            err = r(k) - y;
        else
            [xf, ref_filtered] = get_ss_output(xf, f, r(k));
            err = ref_filtered - y;
        end
        

        x(k) = abs(u_set - u_desirable);
        us(k) = u_set;
        ud(k) = u_desirable;
        if k > 1 && abs(u_set - u_desirable) > anti_windup_tol
            [xc, u] = get_ss_output(xc, c, 0);
            aw(k) = 1;
        else
            [xc, u] = get_ss_output(xc, c, err);
            aw(k) = 0;
        end

        if strcmp(controller_type, 'I+P')
            u = u - (kc * y);
        end

        u_delay = [u u_delay(1:ld)];
        u_desirable = u_delay(ld+1);

        u_set = u_desirable;

        if length(saturation) == 2
            if u_desirable <= saturation(1)
                u_set = saturation(1);
            elseif u_desirable >= saturation(2)
                u_set = saturation(2);
            end
        end

        [xp, y_aux] = get_ss_output(xp, p, (u_set + qu(k)));
        y = y_aux + qy(k);

        U(k) = u_set;
        Y(k) = y;
        E(k) = err;
    end
    hax = subplot(4,1,4);
    plot_signal(sim_time, 'Tempo (s)', x, '', 'm', '$|u_{d}(t) - u(t)|$')
    hold on
    plot_signal(sim_time, 'Tempo (s)', ud, '', 'r', '$u_{d}(t)$')
    plot_signal(sim_time, 'Tempo (s)', us, '', 'b', '$u(t)$')
    plot_signal(sim_time, 'Tempo (s)', aw, 'Valor no instante $t$', 'g', 'anti-windup habilitado')
    set(hax, 'ytick', -0.2:0.2:1.4);
    set(hax, 'ylim', [-0.2 1.4]);
end

printf("Loaded successfuly common functions \n");
