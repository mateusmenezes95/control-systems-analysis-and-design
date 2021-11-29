run base_script.m

%a(2,2) = 0,3(1-xeq1²)
a1 = [0 1; -1 0];
a2 = [0 1; -1 -4.5];
b = [0; 1];
c = [1 0];

dt = integration_step_size;
a1_euler = eye(size(a1)) + dt*a1;
a2_euler = eye(size(a2)) + dt*a2;
b_euler = b*dt;

q = eye(2,2);
r = 0.1;
[k1, x1, l1] = lqr(a1, b, q, r)
[k2, x2, l2] = lqr(a2, b, q, r)

x1(1) = 4;
x2(1) = 0;
ut = 0.1*ones(1, sim_time_length);
xt(:,1) = xt0;

return

for k=1:sim_time_length
    if k != sim_time_length
        x1(k+1) = x2(k)*dt + x1(k);
        x2(k+1) = (-x1(k) + 0.3*(1-(x1(k))^2)*x2(k) + ut(k))*dt + x2(k);
        xt(:,k+1) = a1_euler*xt(:,k) + b_euler*ut(k);
    endif
    y(k) = x1(k);
    yt(k) = c*xt(:,k);
endfor

subplot(2,2,1)
plot_signal(sim.time, 'Tempo [s]', xt(1,:), 'Estado x_1(t)', 'b', 'x_1(t)')
subplot(2,2,2)
plot_signal(sim.time, 'Tempo [s]', xt(2,:), 'Estado x_2(t)', 'b', 'x_2(t)')
subplot(2,2,3:4)
plot_signal(sim.time, 'Tempo [s]', yt, 'Saida y(t)', 'b', 'y(t)')
