run graphs_functions_script

function plot_equilibrium_line(y)
    line(xlim(), [y y], 'color', 'k', 'linestyle', '--')
endfunction

function plot_mimo_references (sim_time, figure_num = 1, reference, legend_str)
    global font_size;
    global line_thickness;
    
    figure(figure_num)

    subplot(2,2,1)
    plot_signal(sim_time, 'x', reference(1,:), 'y', '--k', '$r_1(t)$')
    legend(legend_str);
    set(legend,'location', 'northoutside')
    set(legend,'orientation', 'horizontal')
    hold on;

    subplot(2,2,2);
    set(legend, 'autoupdate', 'off')
    plot_signal(sim_time, '', reference(2,:), '', '--k', '$r_2(t)$')
    set(legend, 'autoupdate', 'on')
    hold on

    subplot(2,2,3);
    plot_signal(sim_time, '', reference(3,:), '', '--k', '$u_{eq1}(t)$')
    legend('off');
    hold on

    subplot(2,2,4);
    plot_signal(sim_time, '', reference(4,:), '', '--k', '$u_{eq2}(t)$')
    legend('off');
    hold on
end

function plot_mimo_response_and_control_signals (sim_time, figure_num = 1,
                                                 y, u, line_color, legend_str)
    global font_size;
    global line_thickness;
    
    figure(figure_num)

    subplot(2,2,1);
    plot_signal(sim_time, 'Tempo (s)',
                y(1,:), 'Saida $y_1(t)$',
                line_color,'')
    legend('off')
    hold on

    subplot(2,2,2);
    plot_signal(sim_time, 'Tempo (s)',
                y(2,:), 'Saida $y_2(t)$',
                line_color, '')

    legend(legend_str);
    set(legend,'location', 'northoutside')
    set(legend,'orientation', 'horizontal')
    hold on

    subplot(2,2,3);
    plot_signal(sim_time, 'Tempo (s)',
                u(1,:), 'Sinal de Controle $u_1(t)$',
                line_color, '')
    legend('off');
    hold on

    subplot(2,2,4);
    plot_signal(sim_time, 'Tempo (s)',
                u(2,:), 'Sinal de Controle $u_2(t)$',
                line_color, '')
    legend('off');
    hold on
end
