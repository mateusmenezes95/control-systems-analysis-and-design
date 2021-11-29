run graphs_functions_script

global line_thickness = 0.8;


function plot_mimo_response_and_control_signals (sim_time, figure_num = 1,
                                                 y, u, reference, 
                                                 line_color, plot_ref = true)
    global font_size;
    global line_thickness;
    
    figure(figure_num)

    subplot(2,2,1);
    if plot_ref
        plot_signal(sim_time, '', reference(1,:), '', '--k', '')
    endif
    hold on
    plot_signal(sim_time, 'Tempo (s)',
                y(1,:), 'Saida $y_1(t)$',
                line_color,'')
    legend('off')

    subplot(2,2,2);
    if plot_ref
        plot_signal(sim_time, '', reference(2,:), '', '--k', '')
    endif
    hold on
    plot_signal(sim_time, 'Tempo (s)',
                y(2,:), 'Saida $y_2(t)$',
                line_color, '')
    legend('off')

    subplot(2,2,3);
    if plot_ref
        plot_signal(sim_time, '', reference(3,:), '', '--k', '')
    endif
    hold on
    plot_signal(sim_time, 'Tempo (s)',
                u(1,:), 'Sinal de Controle $u_1(t)$',
                line_color, '')
    legend('off');

    subplot(2,2,4);
    if plot_ref
        plot_signal(sim_time, '', reference(4,:), '', '--k', '')
    endif
    hold on
    plot_signal(sim_time, 'Tempo (s)',
                u(2,:), 'Sinal de Controle $u_2(t)$',
                line_color, '')
    legend('off');
end

function set_mimo_signals_graph_legend(figure_num, x_pos, legend_str)
    set_figure_size(figure_num)

    subplot(2,2,1)
    legend('show')
    set(legend, 'location', 'northoutside')
    set(legend, 'orientation', 'horizontal')
    legend(legend_str)
    legend_pos = get(legend, 'position');
    legend_pos(1) = x_pos;
    set(legend, 'position', legend_pos);
endfunction

function plot_mimo_states_estimation_error(num_samples, err, line_color, figure_num)
    k = 0:(num_samples - 1);

    global font_size;
    global line_thickness;

    set_figure_size(figure_num)

    for i=1:4
        subplot(2,2,i)
        stem(k, err(i,:), 'color', line_color, 'linewidth', (line_thickness - 0.7), 'marker', 'none')
        ylabel(['Erro $e_' num2str(i) '[k]$'], 'fontsize', font_size);
        xlabel('Amostra [k]', 'fontsize', font_size);
        xlim([0 k(end)])
        ylim([1.2*min(err(i,:)) 1.2*max(err(i,:))])
    endfor
endfunction

function plot_mimo_output_estimation_error(sim_time, err, line_color, figure_num)
    set_figure_size(figure_num)

    for i=1:2
        subplot(2,1,i)
        i_str = num2str(i);
        plot_signal(sim_time, 'Tempo (s)',
                    err(i,:), ['Saida $y_' i_str '(t) - \hat{y}_' i_str '(t)$'],
                    line_color,'')
        ylim([1.2*min(err(i,:)) 1.2*max(err(i,:))])
        legend('off')
        hold on
    endfor
endfunction
