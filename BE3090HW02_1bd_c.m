% Parameters
kgen_values = [0.5, 0.3, 0.15]; % New kgen values for traditional, social media, and individual
kloss_values = [0.657, 1.53, 0.153]; % Corresponding kloss values for traditional, social media, and individual

% New parameters for misinformation
kgen_mis = 0.05; % Generation rate for misinformation
kloss_mis = 0.461; % Loss rate for misinformation
t_mis_injection_duration = 1; % Duration of injection in days
t_mis_decay_duration = 10; % Decay duration for misinformation
injection_interval = 7; % Injection every 7 days
overlap_duration = 3; % Overlap for gradual increase in misinformation

% Time range (0 to 6*30 days)
t = linspace(0, 6*30, 1000);

% Create a figure for plotting
figure;

% Initialize total INF arra3
total_inf = zeros(size(t));

% Loop through the kgen and kloss values for traditional, social media, and individual
for i = 1:length(kgen_values)
    kgen = kgen_values(i);
    kloss = kloss_values(i);
    
    % Function to plot for traditional, social media, and individual
    f_t = (kgen)/kloss * (1 - exp(-kloss * t));
    
    % Add to total INF
    total_inf = total_inf + f_t;
    
    % Plotting the function
    plot(t, f_t, 'LineWidth', 2);
    hold on; % Hold on to plot all functions on the same figure
end

% Plotting misinformation
f_t_mis = zeros(size(t)); % Initialize misinformation contribution array

% Loop for misinformation injection
for day = 0:floor(max(t)/injection_interval)
    % Calculate the time range for the injection
    injection_start = day * injection_interval; % Start of injection
    injection_end = injection_start + t_mis_injection_duration; % End of injection
    
    % Apply the injection (removing the 1- part of the original equation)
    f_t_mis(t >= injection_start & t < injection_end) = kgen_mis * (1 - exp(-kloss_mis * (t(t >= injection_start & t < injection_end) - injection_start)));

    % Calculate decay after the injection period (over the next 10 days)
    decay_start = injection_end; % Start of decay
    decay_end = decay_start + t_mis_decay_duration; % End of decay
    decay_range = t(t >= decay_start & t < decay_end);
    
    % Update misinformation contribution during decay
    f_t_mis(t >= decay_start & t < decay_end) = kgen_mis * exp(-kloss_mis * (decay_range - decay_start));
end

% Add misinformation to total INF
total_inf = total_inf + f_t_mis;

% Plot misinformation
plot(t, f_t_mis, 'LineWidth', 2, 'DisplayName', 'INF misinformation');

% Plot total INF
plot(t, total_inf, 'LineWidth', 2, 'Color', 'k', 'LineStyle', '-.', 'DisplayName', 'Total INF');

% Customize the plot
xlabel('Time (days)');
ylabel('Proportion of INF contribution to INF equi');
title('All INF Sources vs Time with Misinformation');
legend({'INF traditional', ...
         'INF social media', ...
         'INF individual', ...
         'INF misinformation', ...
         'Total INF'}, ...
         'Location', 'Best');
grid on;

% Release hold
hold off;
