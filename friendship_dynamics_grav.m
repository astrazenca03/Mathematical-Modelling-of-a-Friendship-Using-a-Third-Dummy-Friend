function chaotic_friendship_detailed_viz()
    clc; clear; close all;

    %% 1. Setup: Randomized Psychological Stiffness
    disp('Initializing Dual-Reality Simulation with Detailed Legends...');
    
    % --- Step 1: Generate Personalities ---
    num_questions = 20;
    R1 = 2 * rand(num_questions, 1) - 1;
    R2 = 2 * rand(num_questions, 1) - 1;
    R3 = 2 * rand(num_questions, 1) - 1;

    % --- Step 2: Calculate Elastic "Bond Strength" ---
    % Multiplier 35 ensures springs compete with the strong gravity
    K12 = mean(abs(R1 - R2)) * 35; 
    K23 = mean(abs(R2 - R3)) * 35;
    K31 = mean(abs(R3 - R1)) * 35;
    
    K_params = [K12, K23, K31];
    fprintf('Bond Strengths: K12=%.1f, K23=%.1f, K31=%.1f\n', K12, K23, K31);

    %% 2. System Physics
    G = 150;      % Social Gravity (Chaotic driver)
    m = 10;       % Mass (Inertia)
    c = 0.02;     % Damping (Low friction)
    L0 = 4;       % Ideal distance

    %% 3. The Butterfly Effect Setup (Epsilon)
    pos1 = [5, 0, 0];
    pos2 = [-2.5, 4.3, 0];
    pos3 = [-2.5, -4.3, 1]; 
    
    vel1 = [0, 2, 1];
    vel2 = [-1.7, -1, 0];
    vel3 = [1.7, -1, -0.5];

    % Reality A
    y0_A = [pos1, pos2, pos3, vel1, vel2, vel3]';
    
    % Reality B (The Perturbed Timeline)
    epsilon = 1e-3; 
    y0_B = y0_A;
    y0_B(1) = y0_B(1) + epsilon; % Tiny shift in Person 1's X position

    %% 4. Run Simulations
    tspan = 0:0.05:50; 
    options = odeset('RelTol', 1e-6, 'AbsTol', 1e-6);
    
    disp('Simulating Reality A...');
    [t, yA] = ode45(@(t,y) dynamics(t, y, m, G, c, K_params, L0), tspan, y0_A, options);
    
    disp('Simulating Reality B...');
    [~, yB] = ode45(@(t,y) dynamics(t, y, m, G, c, K_params, L0), tspan, y0_B, options);

    %% 5. Visualization
    figure('Name', 'Chaotic Friendship: Butterfly Effect Analysis', 'Color', 'w', 'Position', [50, 50, 1600, 600]);
    
    % --- SUBPLOT 1: REALITY A ---
    subplot(1, 3, 1);
    plot_trajectory(yA, 'Reality A: Baseline Timeline');
    
    % --- SUBPLOT 2: REALITY B ---
    subplot(1, 3, 2);
    plot_trajectory(yB, sprintf('Reality B: Perturbed (\\epsilon = %.3f)', epsilon));
    
    % --- SUBPLOT 3: THE DIVERGENCE ---
    subplot(1, 3, 3);
    
    % Calculate Euclidean Distance
    delta = sqrt(sum((yA - yB).^2, 2));
    
    semilogy(t, delta, 'k-', 'LineWidth', 1.5);
    title('System Divergence (Lyapunov Analysis)');
    xlabel('Time (t)');
    ylabel('State Difference ||A - B|| (Log Scale)');
    grid on;
    
    % Add Legend for Divergence
    legend({'Trajectory Divergence'}, 'Location', 'southeast');
    
    % Fit Lyapunov Trend
    half_idx = round(length(t)/2);
    p = polyfit(t(1:half_idx), log(delta(1:half_idx)), 1);
    lambda = p(1);
    
    % Annotate Graph
    text(t(half_idx), delta(half_idx), sprintf('\\lambda = %.3f', lambda), ...
        'FontSize', 14, 'FontWeight', 'bold', 'Color', 'r', 'BackgroundColor', 'w');
    
    if lambda > 0
        subtitle('Conclusion: CHAOTIC (Sensitive Dependence)');
    else
        subtitle('Conclusion: STABLE (Predictable)');
    end

end

%% Helper: Trajectory Plotter with Detailed Legends
function plot_trajectory(y, plot_title)
    % Extract Coordinates
    P1 = y(:, 1:3);
    P2 = y(:, 4:6);
    P3 = y(:, 7:9);
    
    hold on;
    
    % Plot Paths (Capture handles for legend)
    h1 = plot3(P1(:,1), P1(:,2), P1(:,3), 'r-', 'LineWidth', 0.5); 
    h2 = plot3(P2(:,1), P2(:,2), P2(:,3), 'g-', 'LineWidth', 0.5);
    h3 = plot3(P3(:,1), P3(:,2), P3(:,3), 'b-', 'LineWidth', 0.5);
    
    % Plot End Positions (Markers)
    plot3(P1(end,1), P1(end,2), P1(end,3), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 8);
    plot3(P2(end,1), P2(end,2), P2(end,3), 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 8);
    plot3(P3(end,1), P3(end,2), P3(end,3), 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 8);
    
    % Connect the group at the final moment (The triangle)
    line([P1(end,1) P2(end,1)], [P1(end,2) P2(end,2)], [P1(end,3) P2(end,3)], 'Color', 'k', 'LineStyle', '--');
    line([P2(end,1) P3(end,1)], [P2(end,2) P3(end,2)], [P2(end,3) P3(end,3)], 'Color', 'k', 'LineStyle', '--');
    line([P3(end,1) P1(end,1)], [P3(end,2) P1(end,2)], [P3(end,3) P1(end,3)], 'Color', 'k', 'LineStyle', '--');

    grid on; axis equal;
    
    % Detailed Legend
    legend([h1, h2, h3], ...
           {'Person 1 (Red Trajectory)', ...
            'Person 2 (Green Trajectory)', ...
            'Person 3 (Blue Trajectory)'}, ...
           'Location', 'best', 'FontSize', 8);
           
    % Simplified Axes Labels
    xlabel('X'); ylabel('Y'); zlabel('Z');
    
    title(plot_title);
    view(30, 30); 
    
    % Set fixed limits to ensure easy visual comparison between Plot A and B
    xlim([-20 20]); ylim([-20 20]); zlim([-20 20]); 
end

%% Physics Engine
function dydt = dynamics(~, y, m, G, c, K, L0)
    % Unpack State
    r1 = y(1:3); r2 = y(4:6); r3 = y(7:9);
    v1 = y(10:12); v2 = y(13:15); v3 = y(16:18);
    
    k12 = K(1); k23 = K(2); k31 = K(3);
    
    % Force Calculation 
    
    % Pair 1-2
    d12 = r2 - r1; dist12 = norm(d12); dir12 = d12 / dist12;
    F_grav_12 = (G * m^2 / (dist12^2 + 2)) * dir12; 
    F_spr_12  = k12 * (dist12 - L0) * dir12;
    
    % Pair 2-3
    d23 = r3 - r2; dist23 = norm(d23); dir23 = d23 / dist23;
    F_grav_23 = (G * m^2 / (dist23^2 + 2)) * dir23;
    F_spr_23  = k23 * (dist23 - L0) * dir23;
    
    % Pair 3-1
    d31 = r1 - r3; dist31 = norm(d31); dir31 = d31 / dist31;
    F_grav_31 = (G * m^2 / (dist31^2 + 2)) * dir31;
    F_spr_31  = k31 * (dist31 - L0) * dir31;
    
    % Net Forces
    F1 = F_grav_12 - F_grav_31 + F_spr_12 - F_spr_31 - c*v1;
    F2 = -F_grav_12 + F_grav_23 - F_spr_12 + F_spr_23 - c*v2;
    F3 = -F_grav_23 + F_grav_31 - F_spr_23 + F_spr_31 - c*v3;
    
    dydt = [v1; v2; v3; F1/m; F2/m; F3/m];
end