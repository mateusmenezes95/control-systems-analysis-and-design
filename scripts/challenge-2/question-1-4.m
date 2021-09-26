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

axs = zeros(1,num_plants);

for i = 1:num_plants
    g(i, :) = k(i) .* exp(-s * l(i)) ./ (tau(i) * s + 1);
    delta(i, :) = abs((g(i, :) - gn) ./ gn);
    axs(i) = subplot(num_plants,1,i);
    semilogx(w, delta(i, :));
endfor

lm = max(delta);
linkaxes(axs);
ylim([min(lm) 2.5])

% =============================================================================
% Plot Graphs
% =============================================================================

figure("name", "Análise da Incerteza Multiplicativa")
subplot(3,1,1)
semilogx(w, lm)

comp_sensibility = (c .* gn) ./ (1 + c .* gn);
subplot(3,1,2)
semilogx(w, abs(comp_sensibility))

comp_sensibility_x_lm = abs(comp_sensibility) .* lm;
subplot(3,1,3)
semilogx(w, comp_sensibility_x_lm)
