run base_script.m

for i=1:4
    desirable_eigenvalues(i) = e^(-sampling_period/4);
endfor

kd = place(Gz_ss.A, Gz_ss.B, desirable_eigenvalues)

% Initial contidions
xk0 = zeros(size(xeq_continue));

integration_step_ratio = 30
integration_step_size = sampling_period/integration_step_ratio

a_bar = eye(size(Gs_ss.A))+(integration_step_size*Gs_ss.A);
b_bar = integration_step_size*Gs_ss.B;

gs_ss_for_euler = ss(a_bar, b_bar, Gs_ss.c, Gs_ss.d)

for k=1:length(sim.time)
    if (mod(k, integration_step_ratio) == 1)
        t = (k-1)*dt;
        if k == 1
            uk = -kd*(xk0-xeq_discrete)+ueq;
        else
            uk = -kd*(xk-xeq_discrete)+ueq;
        endif
    endif

    if k == 1
        [xk, yk] = get_ss_output(xk0, gs_ss_for_euler, uk);
        Y = yk;
        U = uk;
        X = xk;
    else
        [xk, yk] = get_ss_output(xk, gs_ss_for_euler, uk);
        Y = [Y yk];
        U = [U uk];
        X = [X xk];
    endif
endfor

y1 = Y(1,:);
y2 = Y(2,:);

u1 = U(1,:);
u2 = U(2,:);

subplot(2,2,1)
plot(sim.time, y1)
set(gca, 'xlim', [0 160])
line(xlim(), [2 2], 'linestyle', '--', 'color', 'k')
ylabel('y_1(t)');
xlabel('Tempo (s)');
subplot(2,2,2)
plot(sim.time, y2)
set(gca, 'xlim', [0 160])
line(xlim(), [1 1], 'linestyle', '--', 'color', 'k')
ylabel('y_2(t)');
xlabel('Tempo (s)');

subplot(2,2,3)
plot(sim.time, u1)
set(gca, 'xlim', [0 160])
line(xlim(), [ueq(1) ueq(1)], 'linestyle', '--', 'color', 'k')
ylabel('u_1(t)');
xlabel('Tempo (s)');
subplot(2,2,4)
plot(sim.time, u2)
set(gca, 'xlim', [0 160])
line(xlim(), [ueq(2) ueq(2)], 'linestyle', '--', 'color', 'k')
ylabel('u_2(t)');
xlabel('Tempo (s)');

eig(Gs_ss.A - Gs_ss.B*kd)