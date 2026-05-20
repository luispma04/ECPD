%% P4 Question 6 

load('singleheater_model_n2.mat')

n = size(A,1);
N = 800;              
t = (0:N-1)*Ts;      

% Open-loop input
u = u_ss*ones(1,N);
Du = u - u_ss; % incremental input, 0

% Real plant simulation variables
Dx = zeros(n,N); 
Dy = zeros(1,N); 
y = zeros(1,N); 

% Kalman filter augmented model
Ad = [A B;
      zeros(1,n) 1];

Bd = [B;
      0];

Cd = [C 0];

%noise
e_std = sqrt(e_var); % measurement noise standard deviation
d_true = 1.15; %input disturbance, bias

% Covariances
QE = Ke * e_var * Ke';
RE = e_var;

deltaE = 1; % tuning parameter, try different 

QEd = [QE zeros(n,1);
       zeros(1,n) deltaE];

% Kalman gain
L = dlqe(Ad, eye(n+1), Cd, QEd, RE);

% Initial estimated augmented state
xd_hat = zeros(n+1,N);

% Initial estimation error, approx 5 degrees Celsius
xd_hat(1:n,1) = pinv(C)*5;

% Estimated variables
Dx_hat = zeros(n,N);
d_hat = zeros(1,N);
Dy_hat = zeros(1,N);
y_hat = zeros(1,N);

Dx_hat(:,1) = xd_hat(1:n,1);
d_hat(1) = xd_hat(end,1);
Dy_hat(1) = C*Dx_hat(:,1);
y_hat(1) = y_ss + Dy_hat(1);

% Initial real output
e = e_std*randn;
Dy(1) = C*Dx(:,1) + e;
y(1) = y_ss + Dy(1);

for k = 1:N-1

    e = e_std*randn;

    Dx(:,k+1) = A*Dx(:,k) + B*(Du(k) + d_true) + Ke*e;
    Dy(k+1) = C*Dx(:,k+1) + e;
    y(k+1) = y_ss + Dy(k+1);

    % Kalman prediction
    xd_pred = Ad*xd_hat(:,k) + Bd*Du(k);

    % Kalman correction
    xd_hat(:,k+1) = xd_pred + L*(Dy(k+1) - Cd*xd_pred);

    % Store estimates
    Dx_hat(:,k+1) = xd_hat(1:n,k+1);
    d_hat(k+1) = xd_hat(end,k+1);

    Dy_hat(k+1) = C*Dx_hat(:,k+1);
    y_hat(k+1) = y_ss + Dy_hat(k+1);

end

% Plots
figure

subplot(2,1,1)
plot(t,y,'LineWidth',1.5)
hold on
plot(t,y_hat,'--','LineWidth',1.5)
grid on
xlabel('Time [s]')
ylabel('Temperature [°C]')
legend('Measured y','Estimated y hat','Location','best')
title('Kalman filter output estimation')

subplot(2,1,2)
plot(t,d_hat,'LineWidth',1.5)
hold on
yline(d_true,'--','LineWidth',1.5)
grid on
xlabel('Time [s]')
ylabel('Disturbance')
legend('Estimated d hat','True d','Location','best')
title('Estimated input disturbance')

sgtitle(['Kalman filter, \delta_E = ' num2str(deltaE)])