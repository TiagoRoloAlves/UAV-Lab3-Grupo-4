clear; clc;

% Parameters
m = 0.035;
g = 9.81;
D = diag([9.1785, 9.1785, 10.311]) * 1e-7;

% Time setup
dt = 0.01;
T_total = 3;
t = 0:dt:T_total;
N = length(t);

% Initial state
x0 = [0.5; -0.3; 0.2; 0; 0; 0];

%% LQR gain
A = [zeros(3), eye(3);
     zeros(3), -D/m];
B = [zeros(3); eye(3)];
Q = diag([100, 100, 100, 10, 10, 10]);
R = eye(3)*2;
K = lqr(A, B, Q, R);

%% --- Linear Simulation using ODE45 ---
[t_lin, x_lin] = ode45(@(t,x) (A - B*K)*x, t, x0);
u_lin = -x_lin * K';

% Extract linear position and velocity
p_lin = x_lin(:,1:3);
v_lin = x_lin(:,4:6);

%% --- Nonlinear Simulation with small-angle R ---
p = zeros(3, N);
v = zeros(3, N);
phi = zeros(1, N);
theta = zeros(1, N);
T_thrust = zeros(1, N);
ua_nl = zeros(N, 3);

p(:,1) = x0(1:3);
v(:,1) = x0(4:6);

for k = 1:N-1
    x = [p(:,k); v(:,k)];
    
    % Compute desired acceleration
    ua = -K * x;
    ua_nl(k,:) = ua';  % store for plot

    % Gravity compensation term
    g_term = [0; 0; -g];
    total_acc = ua - g_term;

    % Recover pitch, roll, thrust
    theta(k) = total_acc(1) / g;
    phi(k) = -total_acc(2) / g;
    T_thrust(k) = m * (total_acc(3) + g);

    % Use small-angle R
    Rsa = small_angle_R(theta(k), phi(k));

    dp = Rsa * v(:,k);
    dv = - (1/m)*D*v(:,k) + ua;

    p(:,k+1) = p(:,k) + dp * dt;
    v(:,k+1) = v(:,k) + dv * dt;
end

theta(end) = theta(end-1);
phi(end) = phi(end-1);
T_thrust(end) = T_thrust(end-1);
ua_nl(end,:) = ua_nl(end-1,:);

%% --- Linear Plots (Figure 1) ---
figure(1);
subplot(3,1,1);
plot(t_lin, p_lin, 'LineWidth', 1.5); title('Linear Position'); ylabel('m'); legend('x','y','z'); grid on;

subplot(3,1,2);
plot(t_lin, v_lin, 'LineWidth', 1.5); title('Linear Velocity'); ylabel('m/s'); legend('v_x','v_y','v_z'); grid on;

subplot(3,1,3);
plot(t_lin, u_lin, 'LineWidth', 1.5); title('Linear Control Input u'); xlabel('Time (s)'); ylabel('u'); legend('u_x','u_y','u_z'); grid on;

%% --- Comparison Plots (Figure 2) ---
figure(2);

% Position
subplot(4,1,1);
plot(t_lin, p_lin, '-', 'LineWidth', 1.2); hold on;
plot(t, p', '--', 'LineWidth', 1.2);
title('Position Comparison'); ylabel('m'); legend('x_{lin}','y_{lin}','z_{lin}','x_{nonlin}','y_{nonlin}','z_{nonlin}'); grid on;

% Velocity
subplot(4,1,2);
plot(t_lin, v_lin, '-', 'LineWidth', 1.2); hold on;
plot(t, v', '--', 'LineWidth', 1.2);
title('Velocity Comparison'); ylabel('m/s'); legend('v_{x,lin}','v_{y,lin}','v_{z,lin}','v_{x,nonlin}','v_{y,nonlin}','v_{z,nonlin}'); grid on;

% ua
subplot(4,1,3);
plot(t_lin, -x_lin*K', '-', 'LineWidth', 1.2); hold on;
plot(t, ua_nl, '--', 'LineWidth', 1.2);
title('u_a Comparison'); ylabel('m/s^2'); legend('u_{x,lin}','u_{y,lin}','u_{z,lin}','u_{x,nonlin}','u_{y,nonlin}','u_{z,nonlin}'); grid on;

% Thrust + pitch + roll
subplot(4,1,4);
plot(t, T_thrust, 'k', 'LineWidth', 1.2); hold on;
plot(t, phi, 'r--', 'LineWidth', 1.2);
plot(t, theta, 'b--', 'LineWidth', 1.2);
title('Nonlinear Thrust and Angles'); xlabel('Time (s)');
ylabel('Thrust (N) / Angle (rad)');
legend('Thrust','Roll (\phi)','Pitch (\theta)'); grid on;

%% --- Small angle rotation matrix
function R = small_angle_R(theta, phi)
    R = [1,     theta*phi,    theta;
         0,     1,           -phi;
        -theta, phi,          1];
end



