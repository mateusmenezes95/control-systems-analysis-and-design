%% -*- texinfo -*-
%% Robust control of a mass-damper-spring system.
##
%% Type @code{which MDSSystem} to locate,
##
%% @code{edit MDSSystem} to open and simply
##
%% @code{MDSSystem} to run the example file.

% ===============================================================================
% Robust Control of a Mass-Damper-Spring System     Lukas Reichlin    August 2011
% ===============================================================================
% Reference: Gu, D.W., Petkov, P.Hr. and Konstantinov, M.M.
%            Robust Control Design with Matlab, Springer 2005

printf("Loading commons functions...\n");

pkg load control
pkg load signal

s = tf("s");

default.dt = 0.01;
default.sim_time = 50;

function time = get_sim_time (dt = default.dt, sim_end_time = default.sim_time)
    time = 0:dt:sim_end_time;
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

printf("Loaded successfuly common functions \n");

Control Systems Analysis And Design
