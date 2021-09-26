addpath(genpath("./")) % Add lib path to Octave script file search paths
run common_parameters_scripts.m

clear s;

w = logspace(-2, 4, 1e4);
s = j*w;

% =============================================================================
% Transfer Functions Definions
% =============================================================================

gn = kn .* exp(-s * ln) ./ (taun * s + 1);
g = zeros(num_plants, length(gn));
delta = zeros(num_plants, length(gn));

c = kc * ((s*ti + 1) ./ s*ti);

f = 1;

% =============================================================================
% Main of the script
% =============================================================================

figure("name", "Incerteza Multiplicativa para cada modelo i")

for i = 1:num_plants
    g(i, :) = k(i) .* exp(-s * l(i)) ./ (tau(i) * s + 1);
    delta(i, :) = abs((g(i, :) - gn) ./ gn);

    axs(i) = subplot(num_plants,1,i);

    semilogx(w, delta(i, :));
    grid on
    set(axs(i), 'ylabel', ['$\Delta_{' num2str(i) '}(jw)$'])
endfor

lm = max(delta);
linkaxes(axs,'y');
ylim([min(lm) 2.5]);
set(axs(1:3), 'xticklabel', []);
set(axs(4), 'xlabel', 'Frequência $w$ (rad/s)');

% =============================================================================
% Plot Graphs
% =============================================================================

clear axs

figure("name", "Análise da Incerteza Multiplicativa")
axs(1) = subplot(3,1,1, 'xticklabel', []);
semilogx(w, lm)
ylabel('$\bar{\Delta}(w)$')
grid on

comp_sensibility = (c .* gn) ./ (1 + c .* gn);
axs(2) = subplot(3,1,2);
semilogx(w, abs(comp_sensibility))
ylabel('$|\textit{C}(jw)|$')
grid on

comp_sensibility_x_lm = abs(comp_sensibility) .* lm;
axs(3) = subplot(3,1,3);
semilogx(w, comp_sensibility_x_lm)
ylabel('$|\textit{C}(jw)|\bar{\Delta}(w)$')
grid on

set(axs(1:2), 'xticklabel', []);
set(axs(3), 'xlabel', 'Frequência $w$ (rad/s)');
