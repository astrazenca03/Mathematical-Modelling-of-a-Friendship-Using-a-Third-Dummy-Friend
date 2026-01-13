function friendship_dynamics()
    clc; clear; close all;

    %% 1. The Questionnaire & Stiffness Generation
    % We simulate inputs for 3 people (P1, P2, P3).
    % In a real scenario, you would use 'input()' to get these values.
    % Scale: [-1, 1]. Magnitude = Stiffness (Stubbornness).
    
    disp('Initializing Psychological Questionnaire...');
    num_questions = 20;
    
    % Generating synthetic responses for demonstration (Simulating User Input)
    % A value of 0.9 means "Very Stubborn/Stiff" on that topic.
    % A value of 0.1 means "Very Flexible".
    responses_P1 = 2 * rand(num_questions, 1) - 1; 
    responses_P2 = 2 * rand(num_questions, 1) - 1;
    responses_P3 = 2 * rand(num_questions, 1) - 1;

    % Create High-Dimensional Stiffness Matrices (20x20)
    % K is diagonal because we assume topics are initially independent
    % We take absolute value because Stiffness must be positive.
    % We add 0.1 so no spring is perfectly broken (0 stiffness).
    K1_full = diag(abs(responses_P1) + 0.1); 
    K2_full = diag(abs(responses_P2) + 0.1);
    K3_full = diag(abs(responses_P3) + 0.1);

    fprintf('Stiffness Matrices created for 3 subjects based on %d questions.\n', num_questions);

    %% 2. Dimensionality Reduction (Compression to 3D)
    % We compress the 20 dimensions into 3 fundamental dimensions (X, Y, Z)
    % to visualize the "Star Configuration" in 3D space.
    % Mapping: Q1-Q7 -> X-axis, Q8-Q14 -> Y-axis, Q15-Q20 -> Z-axis
    
    function K_3D = compress_stiffness(K_full)
        diag_vals = diag(K_full);
        k_x = mean(diag_vals(1:7));   % Average stiffness in "Trait X"
        k_y = mean(diag_vals(8:14));  % Average stiffness in "Trait Y"
        k_z = mean(diag_vals(15:20)); % Average stiffness in "Trait Z"
        K_3D = diag([k_x, k_y, k_z]);
    end

    K1 = compress_stiffness(K1_full);
    K2 = compress_stiffness(K2_full);
    K3 = compress_stiffness(K3_full);

    %% 3. Defining the Dynamic System (Star Configuration)
    % Mass M is connected to 3 walls by springs K1, K2, K3.
    % Equation of Motion: M*x'' + C*x' + K_eq*x = 0
    % where K_eq = K1 + K2 + K3 (Parallel springs holding the center)
    
    M = 5 * eye(3);       % Mass matrix (Inertia of the group)
    C = 2 * eye(3);       % Damping matrix (Emotional regulation/friction)
    K_eq = K1 + K2 + K3;  % Equivalent Stiffness of the star center
    
    % Initial Condition: A "Disagreement" pushes the group away from center
    x0 = [1.0; -0.8; 0.5]; % Initial displacement
    v0 = [0; 0; 0];        % Initial velocity

    % Solve ODE using State Space representation
    % State vector z = [position; velocity]
    tspan = 0:0.1:20;
    [t, z] = ode45(@(t,z) system_dynamics(t, z, M, C, K_eq), tspan, [x0; v0]);
    
    positions = z(:, 1:3); % Extract X, Y, Z coordinates over time

    %% 4. Visualization: 5 Time-Step Snapshots
    figure('Name', 'Friendship Dynamics: The Search for Equilibrium', 'Color', 'w');
    indices = floor(linspace(1, length(t), 5)); % Select 5 evenly spaced time steps
    
    for i = 1:5
        idx = indices(i);
        subplot(1, 5, i);
        
        % Plot the Center (Friendship State)
        plot3(positions(1:idx,1), positions(1:idx,2), positions(1:idx,3), 'b-', 'LineWidth', 1.5); hold on;
        plot3(positions(idx,1), positions(idx,2), positions(idx,3), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 8);
        
        % Plot "Walls" (The 3 Individuals anchoring the friendship)
        % We arbitrarily place them in 3D space to visualize the pull
        scatter3(2, 0, 0, 50, 'k', 'filled'); text(2.1,0,0, 'P1');
        scatter3(-1, 1.7, 0, 50, 'k', 'filled'); text(-1.1,1.8,0, 'P2');
        scatter3(-1, -1.7, 0, 50, 'k', 'filled'); text(-1.1,-1.8,0, 'P3');
        
        % Formatting
        grid on; axis([-2 2 -2 2 -1 1]);
        title(sprintf('T = %.1f s', t(idx)));
        xlabel('Trait X'); ylabel('Trait Y'); zlabel('Trait Z');
        view(45, 30);
    end
    
    sgtitle('Dynamic Convergence of the Friendship Algorithm');

    %% 5. Trend Analysis (Mathematical Formulation)
    % Calculate System Energy (Lyapunov Function) over time
    % V(t) = 0.5 * x' * K_eq * x + 0.5 * v' * M * v
    energy = zeros(length(t), 1);
    for k = 1:length(t)
        pos = z(k, 1:3)';
        vel = z(k, 4:6)';
        energy(k) = 0.5 * pos' * K_eq * pos + 0.5 * vel' * M * vel;
    end
    
    % Fit an exponential decay curve: E(t) = E0 * exp(-lambda * t)
    f = fit(t, energy, 'exp1'); 
    
    % Display Formulation
    fprintf('\n------------------------------------------------\n');
    fprintf('TREND ANALYSIS:\n');
    fprintf('The stability of this friendship follows a Lyapunov Decay.\n');
    fprintf('Mathematical Trend: E(t) = %.2f * e^(%.2f * t)\n', f.a, f.b);
    
    if f.b < -0.5
        fprintf('Conclusion: Highly Stable. Conflicts resolve quickly.\n');
    elseif f.b < 0
        fprintf('Conclusion: Stable, but grudges linger (slow decay).\n');
    else
        fprintf('Conclusion: Unstable. The group cannot reach consensus.\n');
    end
    fprintf('------------------------------------------------\n');
    
    % Plot Energy Trend
    figure('Name', 'Energy Trend', 'Color', 'w');
    plot(f, t, energy);
    title('Dissipation of Conflict Energy over Time');
    xlabel('Time (s)'); ylabel('Psychological Potential Energy');
    grid on;

end

%% Auxiliary Function: System Dynamics ODE
function dzdt = system_dynamics(~, z, M, C, K)
    % State vector z contains [x; y; z; vx; vy; vz]
    pos = z(1:3);
    vel = z(4:6);
    
    % Acceleration = M^-1 * (-C*vel - K*pos)
    acc = M \ (-C * vel - K * pos);
    
    dzdt = [vel; acc];
end