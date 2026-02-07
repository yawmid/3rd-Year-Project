
clc; clear; close all;

% -------------------------------------------------------------------------
% 1. SETUP & DIRECTORIES
% -------------------------------------------------------------------------
fprintf("\n[[ Coronary Artery 4D PCA Analysis ]]\n");
fprintf("Started at %s\n", datetime);

warning('OFF', 'MATLAB:table:ModifiedAndSavedVarnames');

% Parameters
npoint_std = 200;       % Standardized number of points along vessel
numgeometries = 30;     % Number of synthetic arteries to generate

% Define Directories (Update this path if needed)
datadir = 'PATH/TO/YOUR/DIRECTORY'

% Output Directory
synthetic_dir = fullfile(fileparts(fileparts(datadir)), 'synthetic_generated_Stenosed');
if ~exist(synthetic_dir, 'dir')
    mkdir(synthetic_dir);
end

% Get Files
files_all = dir(fullfile(datadir, '*.csv'));
idx = contains({files_all.name}, ".csv");
files_csv = files_all(idx);
npatient = numel(files_csv);

if npatient == 0
    error('No CSV files found in directory: %s', datadir);
end

% -------------------------------------------------------------------------
% 2. LOAD DATA & INTERPOLATE
% -------------------------------------------------------------------------
fprintf("\n[[ Loading and Interpolating Data for %d Patients ]]\n", npatient);

ptdata = struct([]);
interpolated_data = zeros(npatient, npoint_std, 4); % 4 columns: X, Y, Z, R

for ipatient = 1:npatient
    filename = fullfile(files_csv(ipatient).folder, files_csv(ipatient).name);
    
    try
        raw_table = readtable(filename);
        data = table2array(raw_table);
    catch
        data = load(filename);
    end
    
    if size(data, 2) < 4, data(:, 4) = 1.5; end 
    
    
    % Interpolate
    xyz_diff = diff(data(:, 1:3));
    dist_cum = [0; cumsum(sqrt(sum(xyz_diff.^2, 2)))];
    [unique_dist, unique_idx] = unique(dist_cum);
    data_unique = data(unique_idx, :);
    
    spline_func = csaps(unique_dist', data_unique(:, 1:4)', 0.85); %Change parameter for smootjing 
    new_dist = linspace(0, unique_dist(end), npoint_std);
    interp_vals = fnval(spline_func, new_dist)'; 
    
    ptdata(ipatient).cxyzr = interp_vals;
    interpolated_data(ipatient, :, :) = interp_vals;
end

% -------------------------------------------------------------------------
% 3. VISUALISE INPUT DATA 
% -------------------------------------------------------------------------
input_plot_dir = fullfile(synthetic_dir, 'input_plots');
if ~exist(input_plot_dir, 'dir'), mkdir(input_plot_dir); end

for i = 1:min(5, npatient)
    data_to_plot = ptdata(i).cxyzr;
    plot_artery_viz(data_to_plot, sprintf('Input Patient: %s', files_csv(i).name), input_plot_dir, sprintf('input_%02d.png', i));
end

% -------------------------------------------------------------------------
% 4. PCA ANALYSIS
% -------------------------------------------------------------------------
fprintf("\n[[ Conducting 4D PCA ]]\n");

cxyzr_mean = zeros(npoint_std, 4);
cxyzr_coeffs = zeros(npoint_std, 4, 4); 
cxyzr_latent = zeros(npoint_std, 4);    

for ipoint = 1:npoint_std
    xsec_data = squeeze(interpolated_data(:, ipoint, :)); 
    [coeff, ~, latent, ~, ~, mu] = pca(xsec_data);
    
    cxyzr_mean(ipoint, :) = mu;
    n_avail = length(latent);
    cxyzr_coeffs(ipoint, 1:4, 1:n_avail) = coeff;
    cxyzr_latent(ipoint, 1:n_avail) = latent;
    
    % Vector Alignment
    if ipoint > 1
        for ipc = 1:n_avail
            v_curr = cxyzr_coeffs(ipoint, :, ipc);
            v_prev = cxyzr_coeffs(ipoint-1, :, ipc);
            if dot(v_prev, v_curr) < 0
                cxyzr_coeffs(ipoint, :, ipc) = -v_curr;
            end
        end
    end
end

% -------------------------------------------------------------------------
% 5. GENERATE SYNTHETIC GEOMETRIES
% -------------------------------------------------------------------------
fprintf("\n[[ Generating %d Synthetic Geometries ]]\n", numgeometries);

for k = 1:numgeometries
    rng('shuffle');
    weights = normrnd(0, 1.6, 1, 4); 
    
    synth_artery = zeros(npoint_std, 4);
    
    for ipoint = 1:npoint_std
        point_recon = cxyzr_mean(ipoint, :);
        deviation = zeros(1, 4);
        
        for ipc = 1:3 
            vec = reshape(cxyzr_coeffs(ipoint, :, ipc), 1, 4); 
            val = sqrt(cxyzr_latent(ipoint, ipc)); 
            deviation = deviation + (weights(ipc) * val * vec);
        end
        synth_artery(ipoint, :) = point_recon + deviation;
    end
    
    % --- SMOOTHING TO REDUCE WEIRD GEOMETRIES

    % Smooth centerline (X,Y,Z) to remove PCA jitter
    synth_artery(:, 1:3) = smoothdata(synth_artery(:, 1:3), 'gaussian', 30);
    % Smooth radius lightly
    synth_artery(:, 4) = smoothdata(synth_artery(:, 4), 'gaussian', 30);
    
    % Ensure Positive Radius
    synth_artery(:, 4) = max(synth_artery(:, 4), 0.3); 

    generated_artery{k} = synth_artery(:,1:4);

    %Save to .mat file
    save('generated_arteries.mat','generated_artery');
    
    % Save CSV
    csv_filename = fullfile(synthetic_dir, sprintf('synth_centerlines_%02d.csv', k));
    T = array2table(synth_artery, 'VariableNames', {'X', 'Y', 'Z', 'Radius'});
    writetable(T, csv_filename);
    
    % Visualize
    plot_artery_viz(synth_artery, sprintf('Synthetic Geometry %d', k), synthetic_dir, sprintf('synth_viz_%02d.png', k));
    fprintf("Saved Synthetic Geometry %d\n", k);
end

fprintf("\n[[ Analysis Completed at %s ]]\n", datetime);

% -------------------------------------------------------------------------
%  Plot Artery To Visually Analyse Atery
% -------------------------------------------------------------------------

function plot_artery_viz(data_matrix, plot_title, save_dir, save_name)
    cords = data_matrix(:, 1:3);
    radii = data_matrix(:, 4);
    npoint = size(cords, 1);
    
    fig = figure('Visible', 'off'); 
    hold on; grid on; axis equal;
    
    % Plot Centerline
    plot3(cords(:,1), cords(:,2), cords(:,3), 'LineWidth', 2, 'Color', 'b');
    
    % Use gradient for smoother tangent calculation 
    tangents = zeros(npoint, 3);
    tangents(:,1) = gradient(cords(:,1));
    tangents(:,2) = gradient(cords(:,2));
    tangents(:,3) = gradient(cords(:,3));
    
    % Normalise tangents
    tangents = tangents ./ vecnorm(tangents, 2, 2);
    
    theta = linspace(0, 2*pi, 30);
    
    for ip = 1:3:npoint
        radius = radii(ip);
        center = cords(ip, :);
        direction = tangents(ip, :);
        
        % Calculate Rotation Frame
        if abs(direction(3)) > abs(direction(1))
            perp = [1, 0, -direction(1)/direction(3)];
        else
            perp = [0, 1, -direction(2)/direction(3)];
        end
        perp = perp / norm(perp);
        ortho = cross(direction, perp);
        
        R_mat = [perp; ortho; direction];
        
        % Circle Points
        circle_x = radius * cos(theta);
        circle_y = radius * sin(theta);
        circle_z = zeros(size(circle_x));
        
        % Rotate and Translate
        pts = [circle_x; circle_y; circle_z];
        pts_trans = (R_mat' * pts)' + center;
        
        plot3(pts_trans(:,1), pts_trans(:,2), pts_trans(:,3), 'Color', [0.6 0.6 0.6], 'LineWidth', 0.5);
    end
    
    title(plot_title);
    xlabel('X'); ylabel('Y'); zlabel('Z');
    view(3);
    
    saveas(fig, fullfile(save_dir, save_name));
    close(fig);

end
