function friendship_dynamics_interactive()
    clc; clear; close all;

    %% 1. The Questionnaire & Stiffness Generation
    disp('================================================');
    disp('      FRIENDSHIP DYNAMICS SIMULATOR (3-Body)    ');
    disp('================================================');
    disp('Instructions:');
    disp('Please enter responses for 20 questions for each person.');
    disp('Scale: -1 (Strong Disagreement) to 1 (Strong Agreement).');
    disp('0 implies neutral/flexible.');
    disp('Format: Enter as a vector, e.g., [0.5, -0.2, 0.1, ...]');
    disp('------------------------------------------------');

    num_questions = 20;
    
    % --- INPUT FOR PERSON 1 ---
    while true
        fprintf('\nPerson 1: Enter your 20 responses:\n');
        try
            r1 = input('>> ');
            if length(r1) == num_questions && all(r1 >= -1) && all(r1 <= 1)
                responses_P1 = r1(:); % Convert to column vector
                break;
            else
                fprintf('Error: Please enter exactly 20 numbers between -1 and 1.\n');
            end
        catch
            fprintf('Invalid input format. Use brackets [ ].\n');
        end
    end

    % --- INPUT FOR PERSON 2 ---
    while true
        fprintf('\nPerson 2: Enter your 20 responses:\n');
        try
            r2 = input('>> ');
            if length(r2) == num_questions && all(r2 >= -1) && all(r2 <= 1)
                responses_P2 = r2(:);
                break;
            else
                fprintf('Error: Please enter exactly 20 numbers between -1 and 1.\n');
            end
        catch
            fprintf('Invalid input format.\n');
        end
    end

    % --- INPUT FOR PERSON 3 ---
    while true
        fprintf('\nPerson 3: Enter your 20 responses:\n');
        try
            r3 = input('>> ');
            if length(r3) == num_questions && all(r3 >= -1) && all(r3 <= 1)
                responses_P3 = r3(:);
                break;
            else
                fprintf('Error: Please enter exactly 20 numbers between -1 and 1.\n');
            end
        catch
            fprintf('Invalid input format.\n');
        end
    end

    % Create High-Dimensional Stiffness Matrices (20x20)
    % Magnitude = Stiffness. We take abs() because spring constant k > 0.
    % We add 0.05 to ensure no singularity (infinite flexibility).
    K1_full = diag(abs(responses_P1) + 0.05); 
    K2_full = diag(abs(responses_P2) + 0.05);
    K3_full = diag(abs(responses_P3) + 0.05);

    fprintf('\nProcessing Stiffness Matrices... Done.\n');

    %% 2. Dimensionality Reduction (Compression to 3D)
    % Compressing 20D psychology into 3D Physics
    function K_3D = compress_stiffness(K_full)
        diag_vals = diag(K_full);
        % Splitting the 20 questions into 3 "Trait Dimensions"
        k_x = mean(diag_vals(1:7));    % Dimension 1 (e.g., Logic)
        k_y = mean(diag_vals(8:14));   % Dimension 2 (e.g., Emotion)
        k_z = mean(diag_vals(15:20));  % Dimension 3 (e.g., Lifestyle)
        K_3D = diag([k_x, k_y, k_z]);
    end

    K1 = compress_stiffness(K1_full);
    K2 = compress_stiffness(K2_full);
    K3 = compress_stiffness(K3_full);

    %% 3. Defining the Dynamic System (Star Configuration)
    % Physics parameters
    M = 5 * eye(3);       % Mass/Inertia
    C = 1.5 * eye(3);     % Damping/Friction (Social skills)
    K_eq = K1 + K2 + K3;  % Combined Stiffness of the group

    % Initial Disturbance (e.g., an argument)
    x0 = [1.5; -1.5; 1.0]; 
    v0 = [0; 0; 0];

    % Solve ODE
    tspan = 0:0.1:25;
    [t, z] = ode45(@(t,z) system_dynamics(t, z, M, C, K_eq), tspan, [x0; v0]);
    
    positions = z(:, 1:3);

    %% 4. Visualization: 5 Time-Step Snapshots
    figure('Name', 'Friendship Dynamics Analysis', 'Color', 'w');
    indices = floor(linspace(1, length(t), 5)); 
    
    for i = 1:5
        idx = indices(i);
        subplot(1, 5, i);
        
        % Plot Trajectory
        plot3(positions(1:idx,1), positions(1:idx,2), positions(1:idx,3), 'b-', 'LineWidth', 1.0); hold on;
        
        % Plot Current State (Red Ball)
        plot3(positions(idx,1), positions(idx,2), positions(idx,3), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 6);
        
        % Plot Anchor Points (Visual representation of P1, P2, P3 pulls)
        scatter3(2, 0, 0, 30, 'k', 'filled'); 
        scatter3(-1, 1.7, 0, 30, 'k', 'filled'); 
        scatter3(-1, -1.7, 0, 30, 'k', 'filled'); 
        
        grid on; axis([-2 2 -2 2 -2 2]);
        title(sprintf('T = %.1fs', t(idx)));
        xlabel('Dim X'); ylabel('Dim Y'); zlabel('Dim Z');
        view(45, 30);
    end
    
    sgtitle('Dynamics of Friendship Convergence');

    %% 5. Trend Analysis & Mathematical Formulation
    % Calculate Energy V(t)
    energy = zeros(length(t), 1);
    for k = 1:length(t)
        pos = z(k, 1:3)';
        vel = z(k, 4:6)';
        energy(k) = 0.5 * pos' * K_eq * pos + 0.5 * vel' * M * vel;
    end
    
    % Fit Exponential Decay
    f = fit(t, energy, 'exp1'); 
    
    fprintf('\n------------------------------------------------\n');
    fprintf('ANALYSIS RESULTS:\n');
    fprintf('Initial Conflict Energy: %.2f Joules (Psychological Units)\n', f.a);
    fprintf('Decay Rate (Resolution Speed): %.2f\n', f.b);
    fprintf('Trend Equation: E(t) = %.2f * exp(%.2f * t)\n', f.a, f.b);
    
    if f.b < -0.3
        disp('>> STATUS: Healthy Dynamic. Rapid stabilization.');
    elseif f.b < 0
        disp('>> STATUS: Rigid Dynamic. Slow to forgive/forget.');
    else
        disp('>> STATUS: Unstable. Oscillation does not decay.');
    end
    fprintf('------------------------------------------------\n');

end

%% Auxiliary: ODE System
function dzdt = system_dynamics(~, z, M, C, K)
    pos = z(1:3);
    vel = z(4:6);
    acc = M \ (-C * vel - K * pos);
    dzdt = [vel; acc];
end