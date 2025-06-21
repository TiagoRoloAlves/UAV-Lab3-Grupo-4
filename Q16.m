clear;

% System parameters
m = 0.035;
D = diag([9.1785, 9.1785, 10.311]) * 1e-7;

% State-space matrices
A = [zeros(3), eye(3);
     zeros(3), -D/m];

B = [zeros(3);
     eye(3)];

% Cost matrices
Q = diag([100, 100, 100, 10, 10, 10]);
R = eye(3)*2;

% Compute LQR gain
K = lqr(A, B, Q, R);

% Simulation parameters
tspan = 0:0.01:10;
x0 = [0.5; -0.3; 0.2; 0; 0; 0];
x_ref = [1; 1; 1; 0; 0; 0];

% ODE45 simulation
[t, x] = ode45(@(t,x) (A - B*K)*(x - x_ref), tspan, x0);

% Compute control input u = -K(x - x_ref)
x_error = x - x_ref';
u = -x_error * K';

% Plotting all in one figure
figure;

% Position subplot
subplot(3,1,1);
plot(t, x(:,1), 'r', 'LineWidth', 1.5); hold on;
plot(t, x(:,2), 'g', 'LineWidth', 1.5);
plot(t, x(:,3), 'b', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Position (m)');
legend('x', 'y', 'z');
title('LQR Closed-Loop Position Response');
grid on;

% Velocity subplot
subplot(3,1,2);
plot(t, x(:,4), 'r--', 'LineWidth', 1.5); hold on;
plot(t, x(:,5), 'g--', 'LineWidth', 1.5);
plot(t, x(:,6), 'b--', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Velocity (m/s)');
legend('v_x', 'v_y', 'v_z');
title('LQR Closed-Loop Velocity Response');
grid on;

% Input subplot
subplot(3,1,3);
plot(t, u(:,1), 'r:', 'LineWidth', 1.5); hold on;
plot(t, u(:,2), 'g:', 'LineWidth', 1.5);
plot(t, u(:,3), 'b:', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Input (u)');
legend('u_x', 'u_y', 'u_z');
title('Control Input');
grid on;