

t = linspace(0, 60, 600);
center1 = 28;      % peak at 28 s
width1 = 5;        % controls width
amp1 = 1;          % amplitude
signal1 = amp1 * exp(-((t - center1).^2) / (2 * width1^2)) + 0.2;
center2 = 32;      % peak at 32 s
width2 = 5;        % same width for similar shape
amp2 = 0.8;        % slightly lower amplitude
signal2 = amp2 * exp(-((t - center2).^2) / (2 * width2^2)) + 0.15;
figure;
plot(t, signal1, 'LineWidth', 2, 'Color', [0 0.5 0]); % dark green
hold on;
plot(t, signal2, 'LineWidth', 2, 'Color', [0.5 0 0.5]); % dark purple
xlabel('Time (s)');
ylabel('\DeltaF/F');
title('Overlapping Symmetric Calcium Transients');
legend({'Cell 1', 'Cell 2'}, 'Location', 'Northeast');
set(gca, 'FontSize', 14, 'Box', 'off');

%%

t = linspace(0, 60, 600);
center1 = 28;      % peak at 28 s
width1 = 5;        % controls width
amp1 = 1;          % amplitude
signal1 = amp1 * exp(-((t - center1).^2) / (2 * width1^2)) + 0.2;
center2 = 36;      % peak at 36 s (later, more to the right)
width2 = 5;        % same width
amp2 = 0.8;        % slightly lower amplitude
signal2 = amp2 * exp(-((t - center2).^2) / (2 * width2^2)) + 0.15;
figure;
plot(t, signal1, 'LineWidth', 2, 'Color', [0 0.5 0]); % dark green
hold on;
plot(t, signal2, 'LineWidth', 2, 'Color', [0.5 0 0.5]); % dark purple
xlabel('Time (s)');
ylabel('\DeltaF/F');
title('Overlapping Symmetric Calcium Transients');
legend({'Cell 1', 'Cell 2'}, 'Location', 'Northeast');
set(gca, 'FontSize', 14, 'Box', 'off');


%%
t = linspace(0, 60, 600);

% Pink signal
center1 = 28; width1 = 5; amp1 = 1;
signal1 = amp1 * exp(-((t - center1).^2) / (2 * width1^2)) + 0.2;

% Bright yellow signal (shifted closer)
center2 = 33; width2 = 5; amp2_original = 0.8;
signal2_raw = amp2_original * exp(-((t - center2).^2) / (2 * width2^2)) + 0.15;

% Rescale yellow amplitude to match pink
max_signal1 = max(signal1);
max_signal2_raw = max(signal2_raw);
scaling_factor = max_signal1 / max_signal2_raw;
signal2 = (signal2_raw - min(signal2_raw)) * scaling_factor + min(signal2_raw);

figure; hold on;

dotted_len = 15;

% Define windows
pink_start = center1 - dotted_len;
pink_end = center1 + dotted_len;
yellow_start = center2 - dotted_len;
yellow_end = center2 + dotted_len;

% Indices for pink dotted and solid lines
idx_pink_dotted_left = t <= pink_start;
idx_pink_dotted_right = t >= pink_end;
idx_pink_solid = t > pink_start & t < pink_end;

% Indices for yellow dotted and solid lines
idx_yellow_dotted_left = t <= yellow_start;
idx_yellow_dotted_right = t >= yellow_end;
idx_yellow_solid = t > yellow_start & t < yellow_end;

% Overlap interval
overlap_start = max(pink_start, yellow_start);
overlap_end = min(pink_end, yellow_end);
idx_overlap = t >= overlap_start & t <= overlap_end;

% Construct polygon for filled area between curves in overlap
X_fill = [t(idx_overlap), fliplr(t(idx_overlap))];
Y_fill = [signal1(idx_overlap), fliplr(signal2(idx_overlap))];

% Plot filled area between curves (gray, translucent)
fill(X_fill, Y_fill, [0.8 0.8 0.8], 'FaceAlpha', 0.3, 'EdgeColor', 'none');

% Plot pink segments
plot(t(idx_pink_dotted_left), signal1(idx_pink_dotted_left), ':', 'LineWidth', 1.5, 'Color', [1 0.4 0.6]);
plot(t(idx_pink_dotted_right), signal1(idx_pink_dotted_right), ':', 'LineWidth', 1.5, 'Color', [1 0.4 0.6]);
plot(t(idx_pink_solid), signal1(idx_pink_solid), '-', 'LineWidth', 2, 'Color', [1 0.4 0.6]);

% Plot bright yellow segments
plot(t(idx_yellow_dotted_left), signal2(idx_yellow_dotted_left), ':', 'LineWidth', 1.5, 'Color', [1 0.85 0]);
plot(t(idx_yellow_dotted_right), signal2(idx_yellow_dotted_right), ':', 'LineWidth', 1.5, 'Color', [1 0.85 0]);
plot(t(idx_yellow_solid), signal2(idx_yellow_solid), '-', 'LineWidth', 2, 'Color', [1 0.85 0]);

xlabel('Time (s)');
ylabel('\DeltaF/F');
title('Calcium Transients with Shaded Overlap Between Curves');
set(gca, 'FontSize', 14, 'Box', 'off');


%%
t = linspace(0, 60, 600);

% Updated RED: shifted further right for more separation
center1 = 45;  % shifted from 28 -> 36 for further separation
width1 = 5; amp1 = 1;
signal1 = amp1 * exp(-((t - center1).^2) / (2 * width1^2)) + 0.2;

% DARK BLUE: same as before
center2 = 33; width2 = 5; amp2 = 0.8;
signal2_raw = amp2 * exp(-((t - center2).^2) / (2 * width2^2)) + 0.15;

% Rescale blue amplitude to match red
scaling_factor = max(signal1) / max(signal2_raw);
signal2 = (signal2_raw - min(signal2_raw)) * scaling_factor + min(signal2_raw);

figure; hold on;

% Event window ±15s for solid lines
dotted_len = 15;
red_start = center1 - dotted_len; red_end = center1 + dotted_len;
blue_start = center2 - dotted_len; blue_end = center2 + dotted_len;

% Indexing
idx_red_dotted_left = t <= red_start;
idx_red_dotted_right = t >= red_end;
idx_red_solid = t > red_start & t < red_end;

idx_blue_dotted_left = t <= blue_start;
idx_blue_dotted_right = t >= blue_end;
idx_blue_solid = t > blue_start & t < blue_end;

% Overlap region for gray rectangle
overlap_start = max(red_start, blue_start);
overlap_end = min(red_end, blue_end);

% Determine Y limits after dummy plot to ensure correct patch scaling
plot(t, signal1, 'Color', 'none');
plot(t, signal2, 'Color', 'none');
yl = ylim; if yl(1) > 0, yl(1) = 0; end

% Draw gray rectangle
patch([overlap_start overlap_end overlap_end overlap_start], ...
      [yl(1) yl(1) yl(2) yl(2)], ...
      [0.8 0.8 0.8], 'FaceAlpha', 0.3, 'EdgeColor', 'none');

% Plot DARK BLUE first (so it stays in the back)
plot(t(idx_blue_dotted_left), signal2(idx_blue_dotted_left), ':', 'LineWidth', 1.5, 'Color', [0 0 0.5]);
plot(t(idx_blue_dotted_right), signal2(idx_blue_dotted_right), ':', 'LineWidth', 1.5, 'Color', [0 0 0.5]);
plot(t(idx_blue_solid), signal2(idx_blue_solid), '-', 'LineWidth', 2, 'Color', [0 0 0.5]);

% Plot BRIGHT RED second (so it stays in front)
plot(t(idx_red_dotted_left), signal1(idx_red_dotted_left), ':', 'LineWidth', 1.5, 'Color', [1 0 0]);
plot(t(idx_red_dotted_right), signal1(idx_red_dotted_right), ':', 'LineWidth', 1.5, 'Color', [1 0 0]);
plot(t(idx_red_solid), signal1(idx_red_solid), '-', 'LineWidth', 2, 'Color', [1 0 0]);

%%
t = linspace(0, 60, 600);

% BRIGHT RED: shifted further right and lower
center1 = 30; %40        % shifted further right
width1 = 2; amp1 = 0.7;  % lower amplitude to be below blue
signal1 = amp1 * exp(-((t - center1).^2) / (2 * width1^2)) + 0.1; % baseline lower

% DARK BLUE: reference signal
center2 = 25; width2 = 5; amp2 = 1;
signal2 = amp2 * exp(-((t - center2).^2) / (2 * width2^2)) + 0.2;

figure; hold on;

% Event window for solid lines ±15 s
dotted_len = 15;
red_start = 25; %center1 - dotted_len; red_end = center1 + dotted_len; %25
blue_start = 10; %center2 - dotted_len; blue_end = center2 + dotted_len; %10

%Indexing
idx_red_dotted_left = t <= red_start;
idx_red_dotted_right = t >= red_end;
idx_red_solid = t > red_start & t < red_end;

idx_blue_dotted_left = t <= blue_start;
idx_blue_dotted_right = t >= blue_end;
idx_blue_solid = t > blue_start & t < blue_end;

red_end = 35;
% Overlap region
overlap_start = max(red_start, blue_start);
overlap_end = min(red_end, blue_end);

% Determine Y limits dynamically
plot(t, signal1, 'Color', 'none');
plot(t, signal2, 'Color', 'none');
yl = ylim; if yl(1) > 0, yl(1) = 0; end

% Gray rectangle across Y in overlap region
patch([overlap_start overlap_end overlap_end overlap_start], ...
      [yl(1) yl(1) yl(2) yl(2)], ...
      [0.8 0.8 0.8], 'FaceAlpha', 0.3, 'EdgeColor', 'none');

% Plot DARK BLUE first so it stays on top
plot(t(idx_blue_dotted_left), signal2(idx_blue_dotted_left), ':', 'LineWidth', 1.5, 'Color', [0 0 0.5]);
plot(t(idx_blue_dotted_right), signal2(idx_blue_dotted_right), ':', 'LineWidth', 1.5, 'Color', [0 0 0.5]);
plot(t(idx_blue_solid), signal2(idx_blue_solid), '-', 'LineWidth', 2, 'Color', [0 0 0.5]);

% Plot BRIGHT RED second, shifted right and below
plot(t(idx_red_dotted_left), signal1(idx_red_dotted_left), ':', 'LineWidth', 1.5, 'Color', [1 0 0]);
plot(t(idx_red_dotted_right), signal1(idx_red_dotted_right), ':', 'LineWidth', 1.5, 'Color', [1 0 0]);
plot(t(idx_red_solid), signal1(idx_red_solid), '-', 'LineWidth', 2, 'Color', [1 0 0]);

%%
% Time vector
t = linspace(0, 60, 600);

% Define signals
% BRIGHT RED: shifted right and lower
center_red = 30;              % center time
width_red = 2; amp_red = 0.7; % lower amplitude
baseline_red = 0.1;
signal_red = amp_red * exp(-((t - center_red).^2) / (2 * width_red^2)) + baseline_red;

% DARK BLUE: reference signal
center_blue = 25; width_blue = 5; amp_blue = 1;
baseline_blue = 0.2;
signal_blue = amp_blue * exp(-((t - center_blue).^2) / (2 * width_blue^2)) + baseline_blue;

% Event window for solid lines
dotted_len = 15; % window half-length for solid region

% RED event window
red_start = 25;           % fixed start for solid region
red_end = 35;             % fixed end for solid region
idx_red_solid = t > red_start & t < red_end;
idx_red_dotted_left = t <= red_start;
idx_red_dotted_right = t >= red_end;

% BLUE event window
blue_start = 10;          % fixed start for solid region
blue_end = 40;            % fixed end for solid region
idx_blue_solid = t > blue_start & t < blue_end;
idx_blue_dotted_left = t <= blue_start;
idx_blue_dotted_right = t >= blue_end;

% Overlap region for gray rectangle
overlap_start = max(red_start, blue_start);
overlap_end = min(red_end, blue_end);

% Initialize figure
figure; hold on;

% Determine Y limits dynamically
plot(t, signal_red, 'Color', 'none');
plot(t, signal_blue, 'Color', 'none');
yl = ylim; if yl(1) > 0, yl(1) = 0; end

% Plot gray rectangle across Y in overlap region
patch([overlap_start overlap_end overlap_end overlap_start], ...
      [yl(1) yl(1) 1.4 1.4], ...
      [0.8 0.8 0.8], 'FaceAlpha', 0.3, 'EdgeColor', 'none');

% Plot DARK BLUE first (to remain visually behind red)
plot(t(idx_blue_dotted_left), signal_blue(idx_blue_dotted_left), ':', 'LineWidth', 1.5, 'Color', [0 0 0.5]);
plot(t(idx_blue_dotted_right), signal_blue(idx_blue_dotted_right), ':', 'LineWidth', 1.5, 'Color', [0 0 0.5]);
plot(t(idx_blue_solid), signal_blue(idx_blue_solid), '-', 'LineWidth', 2, 'Color', [0 0 0.5]);

% Plot BRIGHT RED second (to remain in front, shifted right and lower)
plot(t(idx_red_dotted_left), signal_red(idx_red_dotted_left), ':', 'LineWidth', 1.5, 'Color', [1 0 0]);
plot(t(idx_red_dotted_right), signal_red(idx_red_dotted_right), ':', 'LineWidth', 1.5, 'Color', [1 0 0]);
plot(t(idx_red_solid), signal_red(idx_red_solid), '-', 'LineWidth', 2, 'Color', [1 0 0]);

%%
% Time vector
t = linspace(0, 60, 600);

% Define signals
% BRIGHT RED: shifted right and lower
center_red = 30; width_red = 2; amp_red = 0.7; baseline_red = 0.1;
signal_red = amp_red * exp(-((t - center_red).^2) / (2 * width_red^2)) + baseline_red;

% DARK BLUE: reference signal
center_blue = 25; width_blue = 5; amp_blue = 1; baseline_blue = 0.2;
signal_blue = amp_blue * exp(-((t - center_blue).^2) / (2 * width_blue^2)) + baseline_blue;

% Event window for solid lines
dotted_len = 15;

% RED event window
red_start = 25;
red_end = 35;
idx_red_solid = t > red_start & t < red_end;
idx_red_dotted_left = t <= red_start;
idx_red_dotted_right = t >= red_end;

% BLUE event window
blue_start = 10;
blue_end = 40;
idx_blue_solid = t > blue_start & t < blue_end;
idx_blue_dotted_left = t <= blue_start;
idx_blue_dotted_right = t >= blue_end;

% Overlap region for gray rectangle
overlap_start = max(red_start, blue_start);
overlap_end = min(red_end, blue_end);

% Initialize figure
figure; hold on;

% Determine Y limits dynamically
plot(t, signal_red, 'Color', 'none');
plot(t, signal_blue, 'Color', 'none');
yl = ylim; if yl(1) > 0, yl(1) = 0; end

% Plot gray rectangle across Y in overlap region
patch([overlap_start overlap_end overlap_end overlap_start], ...
      [yl(1) yl(1) 1.4 1.4], ...
      [0.8 0.8 0.8], 'FaceAlpha', 0.3, 'EdgeColor', 'none');

% Plot DARK BLUE first
plot(t(idx_blue_dotted_left), signal_blue(idx_blue_dotted_left), ':', 'LineWidth', 1.5, 'Color', [0 0 0.5]);
plot(t(idx_blue_dotted_right), signal_blue(idx_blue_dotted_right), ':', 'LineWidth', 1.5, 'Color', [0 0 0.5]);
plot(t(idx_blue_solid), signal_blue(idx_blue_solid), '-', 'LineWidth', 2, 'Color', [0 0 0.5]);

% Plot BRIGHT RED second
plot(t(idx_red_dotted_left), signal_red(idx_red_dotted_left), ':', 'LineWidth', 1.5, 'Color', [1 0 0]);
plot(t(idx_red_dotted_right), signal_red(idx_red_dotted_right), ':', 'LineWidth', 1.5, 'Color', [1 0 0]);
plot(t(idx_red_solid), signal_red(idx_red_solid), '-', 'LineWidth', 2, 'Color', [1 0 0]);

% Annotations for delay and durations
y_arrow = 1.2; % consistent Y position for arrows

% Delay between RED and BLUE onset
delay = red_start - blue_start;
annotation('doublearrow', ...
    [(blue_start/60) (red_start/60)], ...
    [(y_arrow/1.4) (y_arrow/1.4)], ...
    'Color', [0.3 0.3 0.3], 'LineWidth', 1.5);
text(mean([blue_start, red_start]), y_arrow+0.05, ...
    sprintf('Delay = %d s', delay), ...
    'HorizontalAlignment', 'center', 'FontSize', 12);

% Duration of BLUE
annotation('doublearrow', ...
    [(blue_start/60) (blue_end/60)], ...
    [(y_arrow-0.15)/1.4 (y_arrow-0.15)/1.4], ...
    'Color', [0 0 0.5], 'LineWidth', 1.5);
text(mean([blue_start, blue_end]), y_arrow-0.1, ...
    sprintf('%.0f s', blue_end - blue_start), ...
    'HorizontalAlignment', 'center', 'FontSize', 12, 'Color', [0 0 0.5]);

% Duration of RED
annotation('doublearrow', ...
    [(red_start/60) (red_end/60)], ...
    [(y_arrow-0.3)/1.4 (y_arrow-0.3)/1.4], ...
    'Color', [1 0 0], 'LineWidth', 1.5);
text(mean([red_start, red_end]), y_arrow-0.25, ...
    sprintf('%.0f s', red_end - red_start), ...
    'HorizontalAlignment', 'center', 'FontSize', 12, 'Color', [1 0 0]);

%%
% Time vector
t = linspace(0, 60, 600);

% Define signals
% BRIGHT RED: shifted right and lower
center_red = 30; width_red = 2; amp_red = 1; baseline_red = 0.1;
signal_red = amp_red * exp(-((t - center_red).^2) / (2 * width_red^2)) + baseline_red;

% DARK BLUE: reference signal
center_blue = 25; width_blue = 5; amp_blue = 1; baseline_blue = 0.2;
signal_blue = amp_blue * exp(-((t - center_blue).^2) / (2 * width_blue^2)) + baseline_blue;

% Event windows
red_start = 24; red_end = 36;
blue_start = 10; blue_end = 40;

idx_red_solid = t > red_start & t < red_end;
idx_red_dotted_left = t <= red_start;
idx_red_dotted_right = t >= red_end;

idx_blue_solid = t > blue_start & t < blue_end;
idx_blue_dotted_left = t <= blue_start;
idx_blue_dotted_right = t >= blue_end;

% Overlap region for gray rectangle
overlap_start = max(red_start, blue_start);
overlap_end = min(red_end, blue_end);

% Initialize figure
figure; hold on;

% Determine Y limits dynamically
plot(t, signal_red, 'Color', 'none');
plot(t, signal_blue, 'Color', 'none');
yl = ylim; if yl(1) > 0, yl(1) = 0; end

% Plot gray rectangle across Y in overlap region
patch([overlap_start overlap_end overlap_end overlap_start], ...
      [yl(1) yl(1) 1.4 1.4], ...
      [0.8 0.8 0.8], 'FaceAlpha', 0.3, 'EdgeColor', 'none');

% Plot curves
% DARK BLUE first
plot(t(idx_blue_dotted_left), signal_blue(idx_blue_dotted_left), ':', 'LineWidth', 1.5, 'Color', [0 0 0.5]);
plot(t(idx_blue_dotted_right), signal_blue(idx_blue_dotted_right), ':', 'LineWidth', 1.5, 'Color', [0 0 0.5]);
plot(t(idx_blue_solid), signal_blue(idx_blue_solid), '-', 'LineWidth', 2, 'Color', [0 0 0.5]);

% BRIGHT RED second
plot(t(idx_red_dotted_left), signal_red(idx_red_dotted_left), ':', 'LineWidth', 1.5, 'Color', [1 0 0]);
plot(t(idx_red_dotted_right), signal_red(idx_red_dotted_right), ':', 'LineWidth', 1.5, 'Color', [1 0 0]);
plot(t(idx_red_solid), signal_red(idx_red_solid), '-', 'LineWidth', 2, 'Color', [1 0 0]);

% Annotations
ax = gca;
xlim_current = ax.XLim;
ylim_current = ax.YLim;

% Normalized X positions
x_norm = @(x) (x - xlim_current(1)) / (xlim_current(2) - xlim_current(1));
% Normalized Y positions
y_norm = @(y) (y - ylim_current(1)) / (ylim_current(2) - ylim_current(1));

%% Delay arrow (gray) at y = 0.15
delay = red_start - blue_start;
annotation('doublearrow', ...
    [blue_start, red_start], ...
    [0.2,0.2], ...
    'Color', [0.3 0.3 0.3], 'LineWidth', 1.5);
text(mean([blue_start, red_start]), 0.17, ...
    sprintf('Delay'), ...
    'HorizontalAlignment', 'center', 'FontSize', 12);
delay = red_start - blue_start; % should be 15

annotation('doublearrow', ...
    [x_norm(blue_start), x_norm(red_start)], ...
    [y_norm(0.2), y_norm(0.2)], ...
    'Color', [0.3 0.3 0.3], 'LineWidth', 1.5);

text(mean([blue_start, red_start]), 0.17, ...
    'Delay', ...
    'HorizontalAlignment', 'center', 'FontSize', 12);


% Compute normalized units conversion:
ax = gca;
x_limits = xlim(ax);
y_limits = ylim(ax);

x_norm = @(x) (x - x_limits(1)) / (x_limits(2) - x_limits(1));
y_norm = @(y) (y - y_limits(1)) / (y_limits(2) - y_limits(1));

% Create double arrow from X=10 to X=25 at Y=0.2
annotation('doublearrow', ...
    [x_norm(blue_start), x_norm(red_start)], ...
    [y_norm(0.2), y_norm(0.2)], ...
    'Color', [0.3 0.3 0.3], 'LineWidth', 1.5);

% Add label centered between 10 and 25 at Y=0.17 (data coordinates)
text(mean([blue_start, red_start]), 0.17, ...
    'Delay', ...
    'HorizontalAlignment', 'center', 'FontSize', 12);


% Duration of BLUE (y = 0.05)
annotation('doublearrow', ...
    [x_norm(blue_start) x_norm(blue_end)], ...
    [y_norm(0.05) y_norm(0.05)], ...
    'Color', [0 0 0.5], 'LineWidth', 1.5);
text(mean([blue_start, blue_end]), 0.07, ...
    sprintf('%d s', blue_end - blue_start), ...
    'HorizontalAlignment', 'center', 'FontSize', 12, 'Color', [0 0 0.5]);

% Duration of RED (y = 0.15, below red curve)
annotation('doublearrow', ...
    [x_norm(red_start) x_norm(red_end)], ...
    [y_norm(0.15) y_norm(0.15)], ...
    'Color', [1 0 0], 'LineWidth', 1.5);
text(mean([red_start, red_end]), 0.17, ...
    sprintf('%d s', red_end - red_start), ...
    'HorizontalAlignment', 'center', 'FontSize', 12, 'Color', [1 0 0]);

% Labels
xlabel('Time (s)');
ylabel('\DeltaF/F');
title('Calcium Transients with Delay and Duration Annotations');
set(gca, 'FontSize', 14, 'Box', 'off');

