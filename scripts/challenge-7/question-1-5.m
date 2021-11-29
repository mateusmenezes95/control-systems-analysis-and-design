clear all, close all, clc
addpath(genpath("../../lib")) % Add lib path to Octave script file search paths

% =============================================================================
% External packages and script
% =============================================================================

run common_functions_script
run graphs_functions_script

pkg load miscellaneous
format short

% =============================================================================
% Simulation parameters
% =============================================================================

dt = 0.01;
end_time = 50;
sim = get_sim_time (dt, end_time);
sim_time_length = length(sim.time);

% =============================================================================
% Linearization and LQR gain computation
% =============================================================================

%a(2,2) = 0,3(1-xeq1²)
a1 = [0 1; -1 0];
a2 = [0 1; -1 -4.5];
b = [0; 1];
c = [1 0];

eigvec1 = eig(a1);
eigvec2 = eig(a2;

q = eye(2,2);
r = 0.1;
[k1, x1, l1] = lqr(a1, b, q, r);
[k2, x2, l2] = lqr(a2, b, q, r);

% =============================================================================
% Plot function
% =============================================================================

function plot_van_der_pol_oscilator_signals(t, x, y, u, fig_num)
    set_figure_size(fig_num)
    subplot(2,2,1)
    plot_signal(t, 'Tempo [s]', x(1,:), '$Estado x_1(t)$', 'b', '')
    ylim([-6 6])
    legend('off')
    subplot(2,2,3)
    plot_signal(t, 'Tempo [s]', x(2,:), 'Estado $x_2(t)$', 'b', '')
    ylim([-6 6])
    legend('off')
    subplot(2,2,2)
    plot_signal(t, 'Tempo [s]', y, 'Saida y(t)', 'b', '')
    ylim([-6 6])
    legend('off')
    subplot(2,2,4)
    plot_signal(t, 'Tempo [s]', u, 'Sinal de controle $u(t)$', 'b', '')
    legend('off')
endfunction

% =============================================================================
% Simulation of Non Linear model by Euler Forward Approximation
% =============================================================================

xeq1 = [1; 0];
xeq2 = [4; 0];
xeq = [xeq1 xeq2];

kgain = [k1; k2];

uteq(1:2) = [1 4];
xt0 = [-5; -5];

xt(:,1) = xt0;

for i=1:3
    for k=1:sim_time_length
        if i < 3
            ut(k) = -kgain(i,:)*(xt(:,k) - xeq(:,i)) + uteq(i);
        else
            ut(k) = 0;
        endif

        if k != sim_time_length
            xt(:,k+1) = [
                        xt(2,k)*dt + xt(1,k);
                        (-xt(1,k) + 0.3*(1-(xt(1,k))^2)*xt(2,k) + ut(k))*dt + xt(2,k)
                        ]; 
        endif

        yt(k) = c*xt(:,k);
    endfor
    plot_van_der_pol_oscilator_signals(sim.time, xt, yt, ut, i)

    if i == 2
        subplot(2,2,3)
        ylim([-6 9])
    endif
endfor
