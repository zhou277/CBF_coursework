% Coursework 3: Control Barrier Functions and Control Lyapunov Functions for Obstacle Avoidance
% EEEN60122 & 40122 Applied Control and Autonomous Systems (2024-25, 2nd Semester)
% Author: [Your Name]
% Date: [Your Date]

%% Clear Environment
clc;
clear;
close all;

%% Simulation Parameters
dt = 0.01;                % Time step (s)
tf = 20;                  % Total simulation time (s)
n_steps = tf/dt;          % Number of simulation steps

%% Robot Initial State
x = 0; 
y = 0; 
theta = 0; 
v = 0;                    % Initial position, orientation, and velocity

%% Goal Position
goal = [6, 6];            % Target location [x, y]

%% Obstacle Definition
% Single circular obstacle positioned to block the straight-line path.
% [x, y, radius] - Modify based on your and your groupmate's ID as per instructions.
obstacle = [1, 2, 1.5];     

d_min = 0.1;              % Safe distance threshold

%% Control Parameters
v_nom = 1;                % Nominal velocity

% Tuning Parameters (students should tune these three parameters)
penalty_p = 10;           % Slack variable penalty in the QP
gamma1 = 2;               % CLF decay rate for orientation
k_1 = 10;                 % Tuning parameter k1 (affects the safety constraint)

% Fixed Parameters (kept constant for this assignment)
gamma2 = 5 * gamma1;      % CLF decay rate for velocity (fixed relative to gamma1)
k_2 = 1;                  % Tuning parameter k2 (fixed)
lambda = 0.3;            % Weight of the distance-related term in CLF1 relative to the orientation error term (fixed)

%% Simulation Loop
trajectory = zeros(n_steps, 4); % Log: [x, y, theta, v] at each time step
flag = 0;                      % Indicator for early goal arrival

for k = 1:n_steps
    %% Compute Barrier Function (CBF) for Obstacle Avoidance
    % Calculate the relative position to the obstacle.
    dx = x - obstacle(1);
    dy = y - obstacle(2);
    % Barrier function ensures the robot remains at a safe distance from the obstacle.
    B = (dx^2 + dy^2) - obstacle(3)^2 - d_min; 
 
    %% Set Up QP Matrices for Decision Variables [omega; a; delta]
    % H: Quadratic cost matrix penalizing control effort and relaxation variable delta.
    H = diag([1, 2, penalty_p]);
    f = [0; 0; 0];  % Linear cost term is set to zero

    %% CBF Constraint in QP (Safety Constraint)
    % The CBF constraint enforces that the barrier function remains non-negative.
    B_dot = 2 * dx * v * cos(theta) + 2 * dy * v * sin(theta);  % Ḃ(x)

    % LgLfB terms: coefficients of control inputs a and omega
    LgB_omega = 2 * dx * v * sin(theta) - 2 * dy * v * cos(theta);
    LgB_a     = 2 * dx * cos(theta) + 2 * dy * sin(theta);
    
    % Final CBF constraint in the form A_cbf * [omega; a; delta] ≤ b_cbf
    A_cbf = [LgB_omega, LgB_a, 0];
    b_cbf = (k_1+k_2) * B_dot + k_1 * k_2 * B;    

    %% CLF Constraint in QP (Goal Convergence Constraint)
    % Compute distances from the goal.
    dgx = x - goal(1);   % x-distance to goal
    dgy = y - goal(2);   % y-distance to goal
    
    desired_theta = atan(dgy/dgx);  
    e_theta = theta - desired_theta;
    
    % CLF for orientation: combines the angular error and a distance penalty.
    V_1 = (theta - atan(dgy/dgx))^2 + lambda * (dgx^2 + dgy^2);
    % CLF for velocity: penalizes the deviation from the nominal velocity.
    V_2 = (v - v_nom)^2; 
    theta_d_dot = (v*sin(theta)*dgx - v*cos(theta)*dgy)/(dgx^2 +dgy^2);
    % The CLF constraints enforce the exponential decrease of V1 and V2.
    % The relaxation variable delta is included to guarantee feasibility.

    A_clf1 = [2 * e_theta, 0,-1];
    A_clf2 = [0, 2 *(v - v_nom), -1];
    b_clf1 = -gamma1 * V_1 + 2*theta_d_dot * e_theta -lambda*(2*dgx*v*cos(theta) + 2*dgy*v*sin(theta));
    b_clf2 = -gamma2 * (v - v_nom)^2;
    
    % combine
    A_clf = [A_clf1; A_clf2];
    b_clf = [b_clf1; b_clf2];    

    %% Combine CBF and CLF Constraints
    A = [A_cbf; A_clf];
    b = [b_cbf; b_clf];

    %% Solve the Quadratic Program (QP)
    % Decision variables: [omega; a; delta]
    U = quadprog(H, f, A, b, [], [], [], [], [], optimoptions('quadprog', 'Display', 'off'));
    if isempty(U)
        U = [0; 0; 0]; % Default to zero control if the QP is infeasible.
        disp('QP infeasible');
    end
    
    % Extract optimal control inputs.
    omega = U(1);
    a = U(2);
    delta = U(3);

    %% Update Robot State Using the Unicycle Model Dynamics
    v = v + a * dt;
    x = x + v * cos(theta) * dt;
    y = y + v * sin(theta) * dt;
    theta = theta + omega * dt;
    
    % Log the current state.
    trajectory(k, :) = [x, y, theta, v];
    
    %% Check if the Robot Has Reached the Goal
    if norm([x, y] - goal) < 0.1   
        flag = 1; 
        traj_length = k - 1; 
        break;
    end
end

%% Plot Simulation Results
figure;
hold on;

if flag == 1 
    % Plot trajectory up to the point where the goal is reached.
    plot(trajectory(1:traj_length, 1), trajectory(1:traj_length, 2), 'b-', 'LineWidth', 2);
else
    % Plot full trajectory if the goal was not reached within the simulation time.
    plot(trajectory(:, 1), trajectory(:, 2), 'b-', 'LineWidth', 2);
end

% Mark the goal position and draw the obstacle.
scatter(goal(1), goal(2), 100, 'g', 'filled');
viscircles(obstacle(1:2), obstacle(3), 'Color', 'r');

% Add labels, legend, and title.
legend('Robot Path', 'Goal');
xlabel('X Position');
ylabel('Y Position');
title('CBF & CLF-Based Obstacle Avoidance');
grid on;
axis equal;