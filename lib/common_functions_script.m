printf("Loading commons functions...\n");

pkg load control
pkg load signal

s = tf("s");

default.dt = 0.01;
default.sim_time = 50;

function sim = get_sim_time (dt = default.dt, end_time = default.sim_time)
    sim.time_step = dt;
    sim.time = 0:dt:end_time;
    if (nargin == 0)
        printf("No arguments passed. Using default values instead\n")
    end
end

function u = get_step_signal (sim_time, start_time, end_time, amplitude)
    sim_end_time = sim_time(length(sim_time));

    if end_time > sim_end_time
        disp(cstrcat(num2str(end_time), " > ", num2str(sim_end_time), ". Then, end_time = ", num2str(sim_end_time)))
        end_time = sim_end_time;
    end

    u = zeros(length(sim_time), 1);
    start_idx = find(sim_time >= start_time)(1);
    end_idx = find(sim_time >= end_time)(1);
    u(start_idx:end_idx) = amplitude;
end

% function u = get_step_signal (sim_time, param = struct)
%     u = get_step_signal (sim_time, param.start_time, param.end_time, param.amplitude)
% end

function input_signal = get_signal (sim_time, start_time, end_time, amplitude)
    input_signal.start_time = start_time;
    input_signal.end_time = end_time;
    input_signal.amplitude = amplitude;
    input_signal.signal = get_step_signal(sim_time, start_time, end_time, amplitude);
end

% Y(s)/R(s)
function y_to_r_tf = get_output_to_reference_tf (G, C, F)
    y_to_r_tf = minreal(F * ((C * G) / (1 + (C * G))));
    y_to_r_tf.inname = 'R(s)';
    y_to_r_tf.outname = 'Y(s)';
end

% Y(s)/Qu(s)
function y_to_qu_tf = get_output_to_input_disturbance_tf (G, C, F)
    y_to_qu_tf = minreal(F * (G / (1 + (C * G))));
    y_to_qu_tf.inname = 'Qu(s)';
    y_to_qu_tf.outname = 'Y(s)';
end

% Y(s)/Qy(s)
function y_to_qy_tf = get_output_to_output_disturbance_tf (G, C, F)
    y_to_qy_tf = minreal(F * (1 / (1 + (C * G))));
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

printf("Loaded successfuly common functions \n");
