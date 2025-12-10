function Task1_Run_PowerFlow()
clc; close all;

disp("===============================================");
disp("   TASK 1 – THREE-PHASE AC POWER FLOW vs LPF   ");
disp("   (Exact Wang-style matrix using your data)   ");
disp("===============================================");

%% 1) Read your IEEE37_News.xlsx (no artificial changes)
xlsx = 'IEEE37_News.xlsx';
[from,to,Z,loads,N] = read_ieee_excel_3ph(xlsx);

fprintf("Loaded %d-bus, %d-branch system from %s\n", ...
    N, numel(from), xlsx);

%% 2) Build 3-phase Y-bus with mutual coupling
Vbase_LL_kV = 12.66;
Sbase_MVA   = 10;
Y3N = build_ybus_3N(N,from,to,Z,Vbase_LL_kV,Sbase_MVA);

%% 3) Extract ORIGINAL loads and convert to per-unit
P0a_kW = loads.PQperPhase_kWkvar{1};
P0b_kW = loads.PQperPhase_kWkvar{2};
P0c_kW = loads.PQperPhase_kWkvar{3};
Q0a_kvar = loads.PQperPhase_kWkvar{4};
Q0b_kvar = loads.PQperPhase_kWkvar{5};
Q0c_kvar = loads.PQperPhase_kWkvar{6};

Sbase_kVA = Sbase_MVA*1000;

Pa = P0a_kW / Sbase_kVA;
Pb = P0b_kW / Sbase_kVA;
Pc = P0c_kW / Sbase_kVA;
Qa = Q0a_kvar / Sbase_kVA;
Qb = Q0b_kvar / Sbase_kVA;
Qc = Q0c_kvar / Sbase_kVA;

%% 4) Full 3-phase AC power flow (Gauss–Seidel)
Vslack = 1.00;     % 1.0 p.u. slack magnitude (can set 1.02 if you like)

[Va_AC,Vb_AC,Vc_AC] = acpf3N_gs(Y3N,Pa,Pb,Pc,Qa,Qb,Qc,1,Vslack);

A_real = abs(Va_AC);
B_real = abs(Vb_AC);
C_real = abs(Vc_AC);

fprintf('AC PF voltage ranges:\n');
fprintf('  Phase A: [%.4f , %.4f]\n', min(A_real), max(A_real));
fprintf('  Phase B: [%.4f , %.4f]\n', min(B_real), max(B_real));
fprintf('  Phase C: [%.4f , %.4f]\n', min(C_real), max(C_real));

%% 5) Linearized PF (Wang 2017 Eq.(9) style)
lambda = 0.15;   % mix between P–V and Q–V sensitivities

[Va_LPF,Vb_LPF,Vc_LPF] = wang2017_lpf_task1( ...
    Y3N,Pa,Pb,Pc,Qa,Qb,Qc,Vslack,lambda);

A_lpf = abs(Va_LPF);
B_lpf = abs(Vb_LPF);
C_lpf = abs(Vc_LPF);

%% 6) Error metrics
errA = abs(A_real - A_lpf);
errB = abs(B_real - B_lpf);
errC = abs(C_real - C_lpf);

fprintf('\nLPF vs AC PF errors (p.u.):\n');
fprintf('  Phase A: max = %.6f, RMSE = %.6f\n', ...
    max(errA), sqrt(mean(errA.^2)));
fprintf('  Phase B: max = %.6f, RMSE = %.6f\n', ...
    max(errB), sqrt(mean(errB.^2)));
fprintf('  Phase C: max = %.6f, RMSE = %.6f\n', ...
    max(errC), sqrt(mean(errC.^2)));
fprintf('  Overall max error = %.6f\n\n', max([errA;errB;errC]));

%% 7) Plot in the same style as the paper
nodes = 1:N;

figure('Color','w','Position',[100 100 1100 600]); hold on; grid on;
% Real curves
plot(nodes,A_real,'-','LineWidth',2,'Color',[0 0.447 0.741], ...
     'DisplayName','Phase A-Real');
plot(nodes,B_real,'-','LineWidth',2,'Color',[0.85 0.325 0.098], ...
     'DisplayName','Phase B-Real');
plot(nodes,C_real,'-','LineWidth',2,'Color',[0.0 0.6 0.0], ...
     'DisplayName','Phase C-Real');
% LPF curves (dashed with squares)
plot(nodes,A_lpf,'--s','LineWidth',1.8,'MarkerSize',5, ...
     'Color',[0.3 0.7 1.0], 'DisplayName','Phase A-LPF');
plot(nodes,B_lpf,'--s','LineWidth',1.8,'MarkerSize',5, ...
     'Color',[1.0 0.4 0.4], 'DisplayName','Phase B-LPF');
plot(nodes,C_lpf,'--s','LineWidth',1.8,'MarkerSize',5, ...
     'Color',[0.4 0.8 0.8], 'DisplayName','Phase C-LPF');

xlabel('Bus Number','FontSize',13,'FontWeight','bold');
ylabel('Voltage Magnitude (p.u.)','FontSize',13,'FontWeight','bold');
title('Task 1: Three-Phase Voltage (Real vs Linearized DistFlow)', ...
      'FontSize',16,'FontWeight','bold');
legend('Location','best','FontSize',12);

xlim([1 N]);
ylim([min([A_real;B_real;C_real])-0.01 , max([A_real;B_real;C_real])+0.02]);

set(gca,'FontSize',12,'LineWidth',1.1);

disp("Task 1 completed and figure generated.");
disp("===============================================");
end
