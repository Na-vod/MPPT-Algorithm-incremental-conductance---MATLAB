clc;
clear;

%% ================= TIME SAMPLES =================
N = 20;                     % 20 real-time samples
t = 1:N;                    % Time index

%% ================= PV MODULE PARAMETERS =================
Voc_ref = 40;               % Open-circuit voltage at STC (V)
Isc_ref = 9;                % Short-circuit current at STC (A)
G_ref   = 1000;             % Reference irradiance (W/m^2)
T_ref   = 25;               % Reference temperature (°C)

Kv = -0.12;                 % Voltage temp coefficient (V/°C)
Ki =  0.005;                % Current temp coefficient (A/°C)

%% ================= REAL-TIME IRRADIANCE & TEMP =================
G = linspace(600, 1000, N);        % Irradiance variation (W/m^2)
T = linspace(25, 40, N);           % Temperature variation (°C)

%% ================= PV VOLTAGE & CURRENT =================
Vpv = zeros(1,N);
Ipv = zeros(1,N);
Ppv = zeros(1,N);

for k = 1:N
    Ipv(k) = Isc_ref * (G(k)/G_ref) * (1 + Ki*(T(k)-T_ref));
    Voc    = Voc_ref + Kv*(T(k)-T_ref);
    Vpv(k) = 0.8 * Voc;             % Approximate MPP voltage
    Ppv(k) = Vpv(k) * Ipv(k);
end

%% ================= MPPT PARAMETERS =================
D = 0.5;                    % Initial duty cycle
step = 0.01;

V_prev = Vpv(1);
I_prev = Ipv(1);

%% ================= BATTERY PARAMETERS =================
Vbat = 24;                  % Battery voltage (V)
Rbat = 2;                   % Battery resistance (Ohm)

%% ================= STORAGE =================
D_hist = zeros(1,N);
Vout   = zeros(1,N);
Iout   = zeros(1,N);
Pout   = zeros(1,N);

D_hist(1) = D;
Vout(1) = (D/(1-D))*Vpv(1);

%% ================= INC CONDUCTANCE MPPT =================
for k = 2:N

    dV = Vpv(k) - V_prev;
    dI = Ipv(k) - I_prev;

    if dV == 0
        if dI > 0
            D = D - step;
        elseif dI < 0
            D = D + step;
        end
    else
        if (dI/dV) > -(Ipv(k)/Vpv(k))
            D = D - step;
        elseif (dI/dV) < -(Ipv(k)/Vpv(k))
            D = D + step;
        end
    end

    % Limit duty cycle
    D = max(0.05, min(0.95, D));
    D_hist(k) = D;

    % Buck–Boost Converter
    Vout(k) = (D/(1-D)) * Vpv(k);

    % Battery charging
    Iout(k) = (Vout(k) - Vbat) / Rbat;
    if Iout(k) < 0
        Iout(k) = 0;
    end

    Pout(k) = Vout(k) * Iout(k);

    V_prev = Vpv(k);
    I_prev = Ipv(k);
end

%% ================= ONE-LINE OUTPUT =================
fprintf('\nS | G(W/m2) | T(C) | Vpv | Ipv | Ppv | D | Vout | Pbat\n');
fprintf('----------------------------------------------------------\n');

for k = 1:N
    fprintf('%2d | %7.1f | %4.1f | %5.2f | %5.2f | %6.2f | %4.2f | %6.2f | %6.2f\n',...
        k, G(k), T(k), Vpv(k), Ipv(k), Ppv(k), D_hist(k), Vout(k), Pout(k));
end

%% ================= PLOTS =================
figure;

subplot(3,1,1);
plot(Ppv,'b-o');
ylabel('PV Power (W)');
grid on;

subplot(3,1,2);
plot(D_hist,'k-o');
ylabel('Duty Cycle');
grid on;

subplot(3,1,3);
plot(Pout,'g-o');
ylabel('Battery Power (W)');
xlabel('Sample');
grid on;
