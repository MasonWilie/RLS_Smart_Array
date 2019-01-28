close all;
clear all;

%% Parameters
d = 0.5; % Spacing between elements
N = 4; % Number of antenna elements
res = 10000; % Number of steps
M = 600; % Loop iterations

lambda = 0.99; % Forgetting factor

n = 0:N-1;

angle_d = 20; % Desired signal angle in degrees
angle_I1 = 70; % First interfernace signal angle in degrees
angle_I2 = -40; % Second interfernace signal angle in degrees

% Converting angles from degrees to radians
angle_d = deg2rad(angle_d);
angle_I1 = deg2rad(angle_I1);
angle_I2 = deg2rad(angle_I2);


%% Desired Signal
Ts = 1e-3; % Sampling frequency
t = 0:res-1; % Samples
A = 1; % Amplitude

S = A * cos(2 * pi * t * Ts); % Desired signal


%% Interferance Signals

AI1 = 5; % Amplitude of the first interferance angle
AI2 = 5; % Amplitude of the second interferance angle

I1 = AI1 * rand(1, res) - 0.5 * AI1; % Interferance signal 1
I2 = AI2 * rand(1, res) - 0.5 * AI2; % Interferance Signal 2


%% Steering Vectors

a_d =  exp(-1j * 2 * pi * d * n' * sin(angle_d)); % Desired direction
a_i1 = exp(-1j * 2 * pi * d * n' * sin(angle_I1)); % Interferance direction
a_i2 = exp(-1j * 2 * pi * d * n' * sin(angle_I2)); % Interferance direction


%% Creating Recieved Signals

SNR = 100; % Signal to noise ratio of random added noise

NOISE = awgn(ones(N, res), SNR) - 1; % Added noise for each element

x_d  = a_d  * S;
x_i1 = a_i1 * I1;
x_i2 = a_i2 * I2;

x = x_d + x_i1 + x_i2 + NOISE;

%% RLS Algorithm

% Finding the covariance matrix of x for starting point
P = inv(cov(x(:, :)'));
P_n = zeros(N); % Next P

e = zeros(1, M); % Errors
g = 0; % Kalman gain
w = zeros(N, 1); % element weights

space_scalar = round(res / (M+1)) - 1; % Evenly spaced samples


for i = 0:M
    k = i * space_scalar + 1; % Sample number
    e(i + 1) =  S(k) - x(:, k)' * w; % Error Calculation
    g = P * x(:, k) / (lambda + x(:, k)' * P * x(:,k)); % Kalman gain
    P_n = (P - g * x(:, k)' * P) / lambda; % Next R_inv
    w = w + e(i + 1) * g; % Updating weights
   	P = P_n;
end

theta = linspace(-pi / 2, pi / 2, res);

AF = w' * exp(-1j * 2 * pi * d * n' * sin(theta));
AF = AF / max(abs(AF));
AF_dB = mag2db(abs(AF));

figure; hold on; % Plotting array factor
plot(rad2deg(theta), AF_dB);
title('Array Factor');
xlabel('Angle (deg)');
ylabel('Magnitude (dB)');
ylim([-40, 0]);
xlim([-90, 90]);
plot([rad2deg(angle_d), rad2deg(angle_d)], [-1000, 1000], '--');
plot([rad2deg(angle_I1), rad2deg(angle_I1)], [-1000, 1000]);
plot([rad2deg(angle_I2), rad2deg(angle_I2)], [-1000, 1000]);
legend('Desired', 'Interference 1', 'Interference 2');
hold off;

figure; hold on; % Plotting error over time
plot(0:M, e);
xlabel('Iteration');
ylabel('Error');
hold off;






%% 2-D Array
clear all;

%% Parameters
d = 0.5; % Spacing between elements
N = 7; % Number of antenna elements
M = 5;
res = 1000; % Number of steps
iter = 800; % Loop iterations

lambda = 0.99; % Forgetting factor

n = 0:N-1;
m = 0:M-1;

theta_d = 20; % Desired signal angle in degrees
theta_I1 = 70; % First interfernace signal angle in degrees

phi_d = 30;
phi_I1 = 60;

% Converting angles from degrees to radians
theta_d = deg2rad(theta_d);
theta_I1 = deg2rad(theta_I1);

phi_d = deg2rad(phi_d);
phi_I1 = deg2rad(phi_I1);




%% Desired Signal
Ts = 1e-3; % Sampling frequency
t = 0:res-1; % Samples
A = 1; % Amplitude

S = A * cos(2 * pi * t * Ts); % Desired signal


%% Interferance Signals

AI1 = 20; % Amplitude of the first interferance angle

I1 = AI1 * rand(1, res) - 0.5 * AI1; % Interferance signal 1



%% Steering Vectors

a_d_b =  exp(-1j * 2 * pi * d * n' * cos(theta_d)) * ...
    exp(-1j * 2 * pi * d * m' * cos(phi_d))'; % Desired direction
a_i1_b = exp(-1j * 2 * pi * d * n' * cos(theta_I1)) * ...
    exp(-1j * 2 * pi * d * m' * cos(phi_I1))'; % Interferance direction

a_d = zeros(N*M, 1);
a_i1 = zeros(N*M, 1);


for i = 0:N-1
    a_d(i * M + 1:(i + 1) * M) = a_d_b(i + 1, :);
    a_i1(i * M + 1:(i + 1) * M) = a_i1_b(i + 1, :);
    
    
end


%% Creating Recieved Signals

SNR = 120; % Signal to noise ratio of random added noise

NOISE = awgn(ones(N * M, res), SNR) - 1; % Added noise for each element

x_d  = a_d  * S;
x_i1 = a_i1 * I1;


x = x_d + x_i1 + NOISE;

%% RLS Algorithm

% Finding the covariance matrix of x for starting point
P = inv(cov(x(:, :)'));
P_n = zeros(N); % Next P

e = zeros(1, iter); % Errors
g = 0; % Kalman gain
w = zeros(N * M, 1); % element weights

space_scalar = round(res / (iter+1)) - 1; % Evenly spaced samples


for i = 0:iter
    k = i * space_scalar + 1; % Sample number
    e(i + 1) =  S(k) - x(:, k)' * w; % Error Calculation
    g = P * x(:, k) / (lambda + x(:, k)' * P * x(:,k)); % Kalman gain
    P_n = (P - g * x(:, k)' * P) / lambda; % Next R_inv
    w1 = w + e(i + 1) * g; % Updating weights
   	P = P_n;
end

theta = linspace(-pi, pi, res);
phi = linspace(-pi, pi, res);

[phi,theta]=meshgrid(phi,theta);
w = zeros(N, M);

for i = 0:N-1
    w(i + 1, :) = w1(i * M + 1:(i + 1) * M);
end

AF = zeros(res, res);
for i = 1:N
    for j = 1:M
        AF = AF + w(i, j) * exp(-1j * 2 * pi * d * (n(i) * cos(theta)...
            + m(j) * cos(phi)));
    end
end

AF = AF ./ max(max(AF));

[X,Y,Z]=sph2cart(theta,phi,abs(AF));
[Xi1, Yi1, Zi1] = sph2cart(theta_I1, phi_I1, linspace(0, 1.5, res));
[Xd, Yd, Zd] = sph2cart(theta_d, phi_d, linspace(0, 1.5, res));


width = 4;

C = sqrt(X.^2 + Y.^2 + Z.^2);
C = C / max(max(C));

figure; hold on;
title('Array Setup');

for i = 0:N-1
   scatter3((0:(M-1)) * d, zeros(1, M), i * d * ones(1, M), 'filled'); 
end
grid on;
xlabel('X spacing in \lambda');
ylabel('Y spacing in \lambda');
zlabel('Z spacing in \lambda');

hold off;


figure; hold on;
title('Array Factor of 2 Dimensional Array');
xlabel('X');
ylabel('Y');
zlabel('Z');
plot3(Xi1, Yi1, Zi1, 'LineWidth', width, 'Color', 'r');
plot3(Xd, Yd, Zd, 'LineWidth', width, 'Color', 'g');
mesh(X,Y,Z, C);
legend('Interference', 'Desired');
colorbar;
hold off;



