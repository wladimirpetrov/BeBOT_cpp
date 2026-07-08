%% ============================================================
%  3D MPC initial-state trajectory plot for 3 scenarios
%  MATLAB version with:
%   1) original trajectory
%   2) optional smoothed trajectory
%   3) optional extracted points from the smoothed trajectory
%   4) user-specified number of extracted points
%   5) optional saving of extracted points to CSV and MAT
%  ============================================================

clear; clc; close all;

%% ------------------------------------------------------------
% User switches
% -------------------------------------------------------------
plot_original_trajectory = true;
plot_smoothed_trajectory = false;
plot_extracted_points    = false;
label_original_points    = false;
label_extracted_points   = false;
save_extracted_points    = false;

%% ------------------------------------------------------------
% Input CSV from your current 3-scenario caller
% -------------------------------------------------------------
csv_file = 'mpc_for_matlab.csv';
T = readtable(csv_file);

% Make sure trajectory is plotted in iteration order
T = sortrows(T, 'iteration');

%% ------------------------------------------------------------
% Required columns from your CURRENT caller
% -------------------------------------------------------------
required_cols = { ...
    'iteration', ...
    'scenario', ...
    'initial_px', 'initial_py', 'initial_pz', ...
    'target_x', 'target_y', 'target_z', ...
    'cyl_x', 'cyl_y', 'cyl_radius', ...
    'sphere_x', 'sphere_y', 'sphere_z', 'sphere_radius'};

for k = 1:numel(required_cols)
    if ~ismember(required_cols{k}, T.Properties.VariableNames)
        error(['Missing required column in CSV: ', required_cols{k}, newline, ...
               'This plotting script expects the CSV from your current caller, ', ...
               'which saves scenario, target, cylinder, and sphere columns.']);
    end
end

%% ------------------------------------------------------------
% Plot bounds
% -------------------------------------------------------------
x_min = -2.0; x_max =  2.0;
y_min =  0.0; y_max =  3.0;
z_min =  0.0; z_max =  2.0;

x_range = x_max - x_min;
y_range = y_max - y_min;
z_range = z_max - z_min;

%% ------------------------------------------------------------
% Cylinder height for visualization only
% -------------------------------------------------------------
cylinder_height = 1.3;
z_bottom = 0.0;

%% ------------------------------------------------------------
% Smoothing / extraction settings
% -------------------------------------------------------------
num_dense_smooth_points = 2000;   % dense smooth curve resolution
num_extract_points      = 100;     % <<< USER CHOOSES NUMBER OF SAVED POINTS

if num_extract_points < 2
    error('num_extract_points must be at least 2.');
end

%% ------------------------------------------------------------
% Create 3D figure
% -------------------------------------------------------------
fig = figure('Color', 'w', 'Position', [100 100 1000 800]);
ax = axes(fig);
hold(ax, 'on');
grid(ax, 'on');
view(ax, -60, 25);

%% ------------------------------------------------------------
% Plot original trajectory
% -------------------------------------------------------------
if plot_original_trajectory
    plot3(ax, ...
        T.initial_px, ...
        T.initial_py, ...
        T.initial_pz, ...
        '-o', ...
        'LineWidth', 2, ...
        'DisplayName', 'MPC initial states');
end

%% ------------------------------------------------------------
% Build smoothed trajectory and extract equally spaced points
% by arc length
% -------------------------------------------------------------
pts = [T.initial_px, T.initial_py, T.initial_pz];

% Arc-length parameter of original points
ds = sqrt(sum(diff(pts, 1, 1).^2, 2));
s = [0; cumsum(ds)];

% Initialize outputs
x_smooth = [];
y_smooth = [];
z_smooth = [];

x_extract = [];
y_extract = [];
z_extract = [];
s_extract = [];
approx_spacing = NaN;
total_length = 0.0;
extracted_pts = [];

if s(end) > 1e-12
    % --------------------------------------------------------
    % Step 1: dense smooth curve
    % Use pchip for smooth shape with less overshoot than spline
    % --------------------------------------------------------
    s_fine = linspace(0, s(end), num_dense_smooth_points);

    x_smooth = interp1(s, pts(:,1), s_fine, 'pchip');
    y_smooth = interp1(s, pts(:,2), s_fine, 'pchip');
    z_smooth = interp1(s, pts(:,3), s_fine, 'pchip');

    % Plot smoothed curve only if requested
    if plot_smoothed_trajectory
        plot3(ax, ...
            x_smooth, ...
            y_smooth, ...
            z_smooth, ...
            '-', ...
            'LineWidth', 2.5, ...
            'DisplayName', 'Smoothed trajectory');
    end

    % --------------------------------------------------------
    % Step 2: arc length along smoothed curve
    % --------------------------------------------------------
    smooth_pts = [x_smooth(:), y_smooth(:), z_smooth(:)];
    ds_smooth = sqrt(sum(diff(smooth_pts, 1, 1).^2, 2));
    s_smooth = [0; cumsum(ds_smooth)];

    total_length = s_smooth(end);

    % --------------------------------------------------------
    % Step 3: choose EXACT number of extracted points
    % equally spaced in arc length
    % --------------------------------------------------------
    s_extract = linspace(0, total_length, num_extract_points).';

    if num_extract_points > 1
        approx_spacing = total_length / (num_extract_points - 1);
    else
        approx_spacing = 0.0;
    end

    % --------------------------------------------------------
    % Step 4: interpolate extracted points
    % --------------------------------------------------------
    x_extract = interp1(s_smooth, x_smooth, s_extract, 'linear');
    y_extract = interp1(s_smooth, y_smooth, s_extract, 'linear');
    z_extract = interp1(s_smooth, z_smooth, s_extract, 'linear');

    extracted_pts = [x_extract, y_extract, z_extract];

    % --------------------------------------------------------
    % Step 5: save extracted points if requested
    % --------------------------------------------------------
    if save_extracted_points
        extracted_table = table( ...
            (1:size(extracted_pts,1)).', ...
            s_extract, ...
            x_extract, ...
            y_extract, ...
            z_extract, ...
            'VariableNames', {'point_id', 'arc_length', 'x', 'y', 'z'});

        output_csv = 'smoothed_trajectory_extracted_points.csv';
        writetable(extracted_table, output_csv);
        fprintf('Saved extracted points to: %s\n', output_csv);

        output_mat = 'smoothed_trajectory_extracted_points.mat';
        save(output_mat, ...
             'extracted_pts', 'x_extract', 'y_extract', 'z_extract', ...
             's_extract', 'num_extract_points', 'approx_spacing', 'total_length');
        fprintf('Saved extracted points to: %s\n', output_mat);
    end

    % --------------------------------------------------------
    % Step 6: plot extracted points only if requested
    % --------------------------------------------------------
    if plot_extracted_points
        scatter3(ax, ...
            x_extract, ...
            y_extract, ...
            z_extract, ...
            36, ...
            'filled', ...
            'DisplayName', 'Extracted points');

        if label_extracted_points
            for i = 1:length(x_extract)
                text(ax, ...
                    x_extract(i), ...
                    y_extract(i), ...
                    z_extract(i) + 0.02, ...
                    num2str(i), ...
                    'FontSize', 7, ...
                    'HorizontalAlignment', 'center');
            end
        end
    end

    fprintf('Number of extracted points: %d\n', num_extract_points);
    fprintf('Total smoothed trajectory length: %.6f m\n', total_length);
    fprintf('Approximate spacing between extracted points: %.6f m\n', approx_spacing);
else
    warning('Trajectory length is zero or too small. Smoothed/extracted points were not generated.');
end

%% ------------------------------------------------------------
% Label each original MPC iteration
% -------------------------------------------------------------
if label_original_points
    for i = 1:height(T)
        text(ax, ...
            T.initial_px(i), ...
            T.initial_py(i), ...
            T.initial_pz(i) + 0.03, ...
            num2str(T.iteration(i)), ...
            'FontSize', 8, ...
            'HorizontalAlignment', 'left', ...
            'VerticalAlignment', 'bottom');
    end
end

%% ------------------------------------------------------------
% Start point
% -------------------------------------------------------------
scatter3(ax, ...
    T.initial_px(1), ...
    T.initial_py(1), ...
    T.initial_pz(1), ...
    80, ...
    's', ...
    'filled', ...
    'DisplayName', 'Start');

%% ------------------------------------------------------------
% Plot unique targets and scenario-specific obstacles
% -------------------------------------------------------------
scenario_ids = unique(T.scenario);
scenario_ids = sort(scenario_ids);

scenario_rows = zeros(numel(scenario_ids), 1);
for i = 1:numel(scenario_ids)
    scenario_rows(i) = find(T.scenario == scenario_ids(i), 1, 'first');
end

Ts = T(scenario_rows, :);

for i = 1:height(Ts)
    row = Ts(i, :);
    scenario_id_int = row.scenario;

    % ---------------------------------------------------------
    % Cylinder obstacle
    % ---------------------------------------------------------
    if row.cyl_radius > 1e-12
        plot_cylinder_matlab( ...
            ax, ...
            row.cyl_x, ...
            row.cyl_y, ...
            row.cyl_radius, ...
            cylinder_height, ...
            z_bottom);

        text(ax, ...
            row.cyl_x, ...
            row.cyl_y, ...
            z_bottom + cylinder_height + 0.05, ...
            sprintf('S%d cylinder', scenario_id_int), ...
            'FontSize', 8, ...
            'HorizontalAlignment', 'center');
    end

    % ---------------------------------------------------------
    % Sphere obstacle
    % ---------------------------------------------------------
    if row.sphere_radius > 1e-12
        plot_sphere_matlab( ...
            ax, ...
            row.sphere_x, ...
            row.sphere_y, ...
            row.sphere_z, ...
            row.sphere_radius);

        text(ax, ...
            row.sphere_x, ...
            row.sphere_y, ...
            row.sphere_z + row.sphere_radius + 0.05, ...
            sprintf('S%d sphere', scenario_id_int), ...
            'FontSize', 8, ...
            'HorizontalAlignment', 'center');
    end
end

%% ------------------------------------------------------------
% Plot targets after collecting duplicates
% This prevents scenario 2 and scenario 3 same target from being
% plotted twice on top of each other.
% -------------------------------------------------------------
target_mat = round([Ts.target_x, Ts.target_y, Ts.target_z], 10);
[unique_targets, ~, ic] = unique(target_mat, 'rows', 'stable');

for i = 1:size(unique_targets, 1)
    tx = unique_targets(i, 1);
    ty = unique_targets(i, 2);
    tz = unique_targets(i, 3);

    scenario_list = Ts.scenario(ic == i);
    scenario_label = strjoin(string(scenario_list.'), '/');

    scatter3(ax, ...
        tx, ty, tz, ...
        160, ...
        '*', ...
        'DisplayName', sprintf('Target scenario %s', scenario_label));

    text(ax, ...
        tx, ty, tz + 0.08, ...
        sprintf('Target %s', scenario_label), ...
        'FontSize', 9, ...
        'HorizontalAlignment', 'center');
end

%% ------------------------------------------------------------
% Fixed axis limits
% -------------------------------------------------------------
xlim(ax, [x_min, x_max]);
ylim(ax, [y_min, y_max]);
zlim(ax, [z_min, z_max]);

%% ------------------------------------------------------------
% Keep x/y/z physical scale consistent
% -------------------------------------------------------------
pbaspect(ax, [x_range, y_range, z_range]);

%% ------------------------------------------------------------
% Labels and formatting
% -------------------------------------------------------------
title(ax, '3D MPC initial-state trajectory with 3 scenarios');
xlabel(ax, 'x [m]');
ylabel(ax, 'y [m]');
zlabel(ax, 'z [m]');

legend(ax, 'Location', 'best');

%% ------------------------------------------------------------
% Save figure
% -------------------------------------------------------------
output_file = 'mpc_initial_xyz_3_scenarios_scaled.png';
exportgraphics(fig, output_file, 'Resolution', 300);
fprintf('Saved figure to: %s\n', output_file);

%% ============================================================
% Local helper functions
% ============================================================

function plot_cylinder_matlab(ax, x0, y0, radius, height, z_bottom)
    if radius <= 1e-12
        return;
    end

    n_theta = 100;
    n_z = 40;

    theta = linspace(0, 2*pi, n_theta);
    z = linspace(z_bottom, z_bottom + height, n_z);

    [theta_grid, z_grid] = meshgrid(theta, z);

    x_grid = x0 + radius * cos(theta_grid);
    y_grid = y0 + radius * sin(theta_grid);

    surf(ax, ...
        x_grid, y_grid, z_grid, ...
        'FaceAlpha', 0.25, ...
        'EdgeColor', 'none');

    x_circle = x0 + radius * cos(theta);
    y_circle = y0 + radius * sin(theta);

    z_bottom_circle = z_bottom * ones(size(theta));
    z_top_circle = (z_bottom + height) * ones(size(theta));

    plot3(ax, x_circle, y_circle, z_bottom_circle, 'LineWidth', 1.5);
    plot3(ax, x_circle, y_circle, z_top_circle, 'LineWidth', 1.5);

    plot3(ax, ...
        [x0, x0], ...
        [y0, y0], ...
        [z_bottom, z_bottom + height], ...
        ':', ...
        'LineWidth', 1.5);
end

function plot_sphere_matlab(ax, x0, y0, z0, radius)
    if radius <= 1e-12
        return;
    end

    n_u = 80;
    n_v = 40;

    u = linspace(0, 2*pi, n_u);
    v = linspace(0, pi, n_v);

    [u_grid, v_grid] = meshgrid(u, v);

    x_grid = x0 + radius * cos(u_grid) .* sin(v_grid);
    y_grid = y0 + radius * sin(u_grid) .* sin(v_grid);
    z_grid = z0 + radius * cos(v_grid);

    surf(ax, ...
        x_grid, y_grid, z_grid, ...
        'FaceAlpha', 0.30, ...
        'EdgeColor', 'none');

    scatter3(ax, x0, y0, z0, 50, 'x');
end