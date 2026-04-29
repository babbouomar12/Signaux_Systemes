%% run_compare_EKF_UKF_VALEURS_PROF.m
% Comparaison EKF vs UKF pour le TP de poursuite mobile par mesures TOA.
% 
%
% IMPORTANT :
%   - Ce fichier est un SCRIPT complet. Il ne demande aucun argument.
%   - Lancez uniquement ce fichier avec le bouton Run de MATLAB.
%   - Ne lancez pas directement les fonctions locales placees en bas.
%
% Objectif :
%   Comparer EKF et UKF sur exactement les memes donnees :
%       1) meme trajectoire vraie Zmat ;
%       2) memes mesures bruitees Y_noisy ;
%       3) memes matrices Q et R ;
%       4) meme initialisation z0 et P0.
%
% Correction importante par rapport a la version precedente :
%   Dans les figures, on trace d'abord l'UKF puis l'EKF au-dessus.
%   Cela evite que la courbe EKF rouge pointillee soit cachee par la courbe UKF bleue.
%
% Figures generees :
%   fig_01_trajectoire_EKF_vs_UKF.png
%   fig_02_EQM_EKF_vs_UKF.png
%   fig_03_DtoP_EKF_vs_UKF.png
%   fig_04_traceP_EKF_vs_UKF.png
%   fig_05_MC_EQM_moyenne_EKF_vs_UKF.png
%   fig_06_MC_DtoP_moyen_EKF_vs_UKF.png
%   fig_07_MC_traceP_moyenne_EKF_vs_UKF.png
%   fig_08_difference_traceP_UKF_moins_EKF.png
%
% Convention :
%   sigmau et sigmav sont des VARIANCES.

clear; close all; clc;

fprintf('============================================================\n');
fprintf('Comparaison EKF vs UKF - Poursuite mobile par TOA\n');
fprintf('Convention : sigmau et sigmav sont des VARIANCES.\n');
fprintf('============================================================\n\n');

%% ========================================================================
%% 1) Dossier de sauvegarde des figures
%% ========================================================================

outDir = fullfile(pwd, 'resultats_EKF_vs_UKF_VALEURS_PROF');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

%% ========================================================================
%% 2) Positions des stations de base
%% ========================================================================

Mk1 = [-1500; 1500];
Mk2 = [20000; 20000];
Mk3 = [25000; 500];

Anchor_position = [Mk1 Mk2 Mk3];
Ns = size(Anchor_position, 2);

%% ========================================================================
%% 3) Parametres de simulation communs
%% ========================================================================

Ts = 1;                  % periode d'echantillonnage
N = 500;                 % nombre d'iterations

sigmau = 2;              % variance du bruit d'acceleration
sigmav = 1e-12;          % variance du bruit de mesure TOA

Zinit = [1000; 1000; 0; 0];
S_seed = 50;

epsilon = 0.01;          % VALEUR PROF C.1 : P(0|0) = epsilon * I4
z0 = [0; 0; 0; 0];
P0 = epsilon * eye(4);

% Valeurs imposees par l'enonce pour la comparaison directe C.1 :
%   generate_traj([1000;1000;0;0], 2, 1, 500, 50, Anchor_position)
%   sigmav = 1e-12, Nw = 1, z_est = [0;0;0;0], P_est = 0.01*I4
% Pour la comparaison Monte-Carlo de type C.3, l'enonce fixe epsilon = 100
% et Nw = 20. On utilisera donc P0_MC = 100*I4 dans la partie MC.

eps_seuil = 100;         % seuil utilise pour detecter la convergence

%% ========================================================================
%% 4) Parametres UKF
%% ========================================================================
% alpha controle la dispersion des points sigma.
% beta = 2 est le choix classique si l'etat est suppose gaussien.
% kappa = 3 est garde pour rester coherent avec l'analyse de l'impact de alpha.

alpha = 0.1;
beta = 2;
kappa = 3;

%% ========================================================================
%% 5) Generation de la trajectoire vraie et des TOA non bruitees
%% ========================================================================

[Zmat, Yf] = generate_traj_local(Zinit, sigmau, Ts, N, S_seed, Anchor_position);

%% ========================================================================
%% 6) Bruitage unique des TOA pour comparaison directe Nw = 1
%% ========================================================================
% Meme Y_noisy pour EKF et UKF : comparaison equitable.

rng(100);
Y_noisy = Yf + sqrt(sigmav) * randn(Ns, N);

%% ========================================================================
%% 7) Lancement EKF et UKF sur les memes donnees
%% ========================================================================

[Z_EKF, EQM_EKF, DtoP_EKF, traceP_EKF] = EKF_track_traj_local( ...
    Zmat, Y_noisy, z0, P0, sigmau, sigmav, Ts, N, Anchor_position);

[Z_UKF, EQM_UKF, DtoP_UKF, traceP_UKF] = UKF_track_traj_local( ...
    Zmat, Y_noisy, z0, P0, sigmau, sigmav, Ts, N, ...
    Anchor_position, alpha, beta, kappa);

%% ========================================================================
%% 8) Diagnostic numerique Nw = 1
%% ========================================================================

fprintf('\n--- Diagnostic Nw = 1 : differences entre EKF et UKF ---\n');
fprintf('max |EQM_UKF    - EQM_EKF|    = %.6e\n', max(abs(EQM_UKF - EQM_EKF)));
fprintf('max |DtoP_UKF   - DtoP_EKF|   = %.6e\n', max(abs(DtoP_UKF - DtoP_EKF)));
fprintf('max |traceP_UKF - traceP_EKF| = %.6e\n', max(abs(traceP_UKF - traceP_EKF)));

%% ========================================================================
%% 9) Figures comparatives Nw = 1
%% ========================================================================

% ------------------------------------------------------------------------
% Figure 1 : Trajectoire
% ------------------------------------------------------------------------
figure('Name','Trajectoire vraie, EKF et UKF');
hold on; grid on;

plot(Zmat(1,:), Zmat(2,:), 'k-', ...
    'LineWidth', 2.2, 'DisplayName', 'Trajectoire vraie');

% On trace UKF avant EKF pour que EKF reste visible au-dessus.
plot(Z_UKF(1,:), Z_UKF(2,:), 'b-.', ...
    'LineWidth', 1.6, 'DisplayName', 'Estimation UKF');

plot(Z_EKF(1,:), Z_EKF(2,:), 'r--', ...
    'LineWidth', 2.2, 'DisplayName', 'Estimation EKF');

plot(Anchor_position(1,:), Anchor_position(2,:), 'ks', ...
    'MarkerFaceColor', 'y', 'MarkerSize', 8, 'DisplayName', 'Stations TOA');

xlabel('x (m)');
ylabel('y (m)');
title('Comparaison de trajectoire : EKF vs UKF');
legend('Location','best');
axis equal;

saveas(gcf, fullfile(outDir, 'fig_01_trajectoire_EKF_vs_UKF.png'));

% ------------------------------------------------------------------------
% Figure 2 : EQM
% ------------------------------------------------------------------------
figure('Name','EQM EKF vs UKF');
hold on; grid on;

plot(EQM_UKF, 'b-', ...
    'LineWidth', 1.6, 'DisplayName', 'EQM UKF');

plot(EQM_EKF, 'r--', ...
    'LineWidth', 2.2, 'DisplayName', 'EQM EKF');

xlabel('Iteration n');
ylabel('EQM');
title('Comparaison de l''Erreur Quadratique Moyenne : EKF vs UKF');
legend('Location','best');

saveas(gcf, fullfile(outDir, 'fig_02_EQM_EKF_vs_UKF.png'));

% ------------------------------------------------------------------------
% Figure 3 : DtoP
% ------------------------------------------------------------------------
figure('Name','DtoP EKF vs UKF');
hold on; grid on;

plot(DtoP_UKF, 'b-', ...
    'LineWidth', 1.6, 'DisplayName', 'DtoP UKF');

plot(DtoP_EKF, 'r--', ...
    'LineWidth', 2.2, 'DisplayName', 'DtoP EKF');

xlabel('Iteration n');
ylabel('DtoP moyen (m)');
title('Comparaison de l''erreur moyenne de position : EKF vs UKF');
legend('Location','best');

saveas(gcf, fullfile(outDir, 'fig_03_DtoP_EKF_vs_UKF.png'));

% ------------------------------------------------------------------------
% Figure 4 : trace(P)
% ------------------------------------------------------------------------
figure('Name','Trace de P EKF vs UKF');
hold on; grid on;

plot(traceP_UKF, 'b-', ...
    'LineWidth', 1.6, 'DisplayName', 'trace(P) UKF');

plot(traceP_EKF, 'r--', ...
    'LineWidth', 2.6, 'DisplayName', 'trace(P) EKF');

xlabel('Iteration n');
ylabel('trace(P)');
title('Comparaison de la covariance estimee : EKF vs UKF');
legend('Location','best');

saveas(gcf, fullfile(outDir, 'fig_04_traceP_EKF_vs_UKF.png'));

%% ========================================================================
%% 10) Tableau numerique pour Nw = 1
%% ========================================================================

[nc_EKF, EQM_final_EKF, DtoP_final_EKF] = final_perf_local(EQM_EKF, DtoP_EKF, eps_seuil);
[nc_UKF, EQM_final_UKF, DtoP_final_UKF] = final_perf_local(EQM_UKF, DtoP_UKF, eps_seuil);

fprintf('\n--- Comparaison directe sur une seule realisation ---\n');
fprintf('Filtre | nc  | EQM_final      | DtoP_final (m) | EQM_fin instant | DtoP_fin instant\n');
fprintf('EKF    | %3d | %13.6e | %14.6f | %15.6e | %16.6f\n', ...
    nc_EKF, EQM_final_EKF, DtoP_final_EKF, EQM_EKF(end), DtoP_EKF(end));
fprintf('UKF    | %3d | %13.6e | %14.6f | %15.6e | %16.6f\n', ...
    nc_UKF, EQM_final_UKF, DtoP_final_UKF, EQM_UKF(end), DtoP_UKF(end));

%% ========================================================================
%% 11) Comparaison Monte-Carlo Nw = 20
%% ========================================================================

Nw = 20;                 % VALEUR PROF C.3
epsilon_MC = 100;          % VALEUR PROF C.3
P0_MC = epsilon_MC * eye(4);

EQM_EKF_mat = zeros(Nw, N);
EQM_UKF_mat = zeros(Nw, N);

DtoP_EKF_mat = zeros(Nw, N);
DtoP_UKF_mat = zeros(Nw, N);

traceP_EKF_mat = zeros(Nw, N);
traceP_UKF_mat = zeros(Nw, N);

for nw = 1:Nw
    rng(1000 + nw);

    % Meme bruit pour EKF et UKF dans cette realisation Monte-Carlo.
    Y_noisy_mc = Yf + sqrt(sigmav) * randn(Ns, N);

    [~, EQM_e, DtoP_e, traceP_e] = EKF_track_traj_local( ...
        Zmat, Y_noisy_mc, z0, P0_MC, sigmau, sigmav, Ts, N, Anchor_position);

    [~, EQM_u, DtoP_u, traceP_u] = UKF_track_traj_local( ...
        Zmat, Y_noisy_mc, z0, P0_MC, sigmau, sigmav, Ts, N, ...
        Anchor_position, alpha, beta, kappa);

    EQM_EKF_mat(nw,:) = EQM_e;
    EQM_UKF_mat(nw,:) = EQM_u;

    DtoP_EKF_mat(nw,:) = DtoP_e;
    DtoP_UKF_mat(nw,:) = DtoP_u;

    traceP_EKF_mat(nw,:) = traceP_e;
    traceP_UKF_mat(nw,:) = traceP_u;
end

EQM_EKF_moy = mean(EQM_EKF_mat, 1);
EQM_UKF_moy = mean(EQM_UKF_mat, 1);

DtoP_EKF_moy = mean(DtoP_EKF_mat, 1);
DtoP_UKF_moy = mean(DtoP_UKF_mat, 1);

traceP_EKF_moy = mean(traceP_EKF_mat, 1);
traceP_UKF_moy = mean(traceP_UKF_mat, 1);

%% ========================================================================
%% 12) Diagnostic numerique Monte-Carlo
%% ========================================================================

diff_EQM_MC = EQM_UKF_moy - EQM_EKF_moy;
diff_DtoP_MC = DtoP_UKF_moy - DtoP_EKF_moy;
diff_traceP_MC = traceP_UKF_moy - traceP_EKF_moy;

fprintf('\n--- Diagnostic Monte-Carlo Nw = %d ---\n', Nw);
fprintf('max |EQM_UKF_moy    - EQM_EKF_moy|    = %.6e\n', max(abs(diff_EQM_MC)));
fprintf('mean |EQM_UKF_moy   - EQM_EKF_moy|    = %.6e\n', mean(abs(diff_EQM_MC)));
fprintf('max |DtoP_UKF_moy   - DtoP_EKF_moy|   = %.6e\n', max(abs(diff_DtoP_MC)));
fprintf('mean |DtoP_UKF_moy  - DtoP_EKF_moy|   = %.6e\n', mean(abs(diff_DtoP_MC)));
fprintf('max |traceP_UKF_moy - traceP_EKF_moy| = %.6e\n', max(abs(diff_traceP_MC)));
fprintf('mean |traceP_UKF_moy- traceP_EKF_moy| = %.6e\n', mean(abs(diff_traceP_MC)));

%% ========================================================================
%% 13) Figures comparatives Monte-Carlo
%% ========================================================================

% ------------------------------------------------------------------------
% Figure 5 : EQM moyenne
% ------------------------------------------------------------------------
figure('Name','Monte-Carlo EQM moyenne EKF vs UKF');
hold on; grid on;

plot(EQM_UKF_moy, 'b-', ...
    'LineWidth', 1.8, 'DisplayName', 'Moyenne UKF');

plot(EQM_EKF_moy, 'r--', ...
    'LineWidth', 2.8, 'DisplayName', 'Moyenne EKF');

xlabel('Iteration n');
ylabel('EQM moyenne');
title(['Monte-Carlo Nw = ' num2str(Nw) ' : EQM moyenne EKF vs UKF']);
legend('Location','best');

saveas(gcf, fullfile(outDir, 'fig_05_MC_EQM_moyenne_EKF_vs_UKF.png'));

% ------------------------------------------------------------------------
% Figure 6 : DtoP moyen
% ------------------------------------------------------------------------
figure('Name','Monte-Carlo DtoP moyen EKF vs UKF');
hold on; grid on;

plot(DtoP_UKF_moy, 'b-', ...
    'LineWidth', 1.8, 'DisplayName', 'Moyenne UKF');

plot(DtoP_EKF_moy, 'r--', ...
    'LineWidth', 2.8, 'DisplayName', 'Moyenne EKF');

xlabel('Iteration n');
ylabel('DtoP moyen (m)');
title(['Monte-Carlo Nw = ' num2str(Nw) ' : DtoP moyen EKF vs UKF']);
legend('Location','best');

saveas(gcf, fullfile(outDir, 'fig_06_MC_DtoP_moyen_EKF_vs_UKF.png'));

% ------------------------------------------------------------------------
% Figure 7 : trace(P) moyenne
% ------------------------------------------------------------------------
figure('Name','Monte-Carlo trace(P) moyenne EKF vs UKF');
hold on; grid on;

% Correction de visibilite :
% On trace UKF d'abord, puis EKF au-dessus.
plot(traceP_UKF_moy, 'b-', ...
    'LineWidth', 1.8, 'DisplayName', 'Moyenne UKF');

plot(traceP_EKF_moy, 'r--', ...
    'LineWidth', 3.0, 'DisplayName', 'Moyenne EKF');

xlabel('Iteration n');
ylabel('trace(P) moyenne');
title(['Monte-Carlo Nw = ' num2str(Nw) ' : trace(P) moyenne EKF vs UKF']);
legend('Location','best');

saveas(gcf, fullfile(outDir, 'fig_07_MC_traceP_moyenne_EKF_vs_UKF.png'));

% ------------------------------------------------------------------------
% Figure 8 : difference trace(P)
% ------------------------------------------------------------------------
figure('Name','Difference trace(P) UKF - EKF');
grid on;

plot(diff_traceP_MC, 'k-', 'LineWidth', 1.8);

xlabel('Iteration n');
ylabel('trace(P) UKF - trace(P) EKF');
title('Difference entre les covariances moyennes UKF et EKF');

saveas(gcf, fullfile(outDir, 'fig_08_difference_traceP_UKF_moins_EKF.png'));

%% ========================================================================
%% 14) Performances finales Monte-Carlo
%% ========================================================================

[nc_EKF_mc, EQM_final_EKF_mc, DtoP_final_EKF_mc] = final_perf_local(EQM_EKF_moy, DtoP_EKF_moy, eps_seuil);
[nc_UKF_mc, EQM_final_UKF_mc, DtoP_final_UKF_mc] = final_perf_local(EQM_UKF_moy, DtoP_UKF_moy, eps_seuil);

fprintf('\n--- Comparaison Monte-Carlo C.3 : epsilon = 100, Nw = %d ---\n', Nw);
fprintf('Filtre | nc  | EQM_final moyen | DtoP_final moyen (m)\n');
fprintf('EKF    | %3d | %15.6e | %20.6f\n', nc_EKF_mc, EQM_final_EKF_mc, DtoP_final_EKF_mc);
fprintf('UKF    | %3d | %15.6e | %20.6f\n', nc_UKF_mc, EQM_final_UKF_mc, DtoP_final_UKF_mc);

%% ========================================================================
%% 15) Sauvegarde d'un resume CSV
%% ========================================================================

csvFile = fullfile(outDir, 'resume_comparaison_EKF_UKF.csv');
fid = fopen(csvFile, 'w');

fprintf(fid, 'Cas,Filtre,nc,EQM_final,DtoP_final_m,EQM_fin_instant,DtoP_fin_instant_m\n');

fprintf(fid, 'Nw1,EKF,%d,%.10e,%.10f,%.10e,%.10f\n', ...
    nc_EKF, EQM_final_EKF, DtoP_final_EKF, EQM_EKF(end), DtoP_EKF(end));

fprintf(fid, 'Nw1,UKF,%d,%.10e,%.10f,%.10e,%.10f\n', ...
    nc_UKF, EQM_final_UKF, DtoP_final_UKF, EQM_UKF(end), DtoP_UKF(end));

fprintf(fid, 'MC%d,EKF,%d,%.10e,%.10f,%.10e,%.10f\n', ...
    Nw, nc_EKF_mc, EQM_final_EKF_mc, DtoP_final_EKF_mc, EQM_EKF_moy(end), DtoP_EKF_moy(end));

fprintf(fid, 'MC%d,UKF,%d,%.10e,%.10f,%.10e,%.10f\n', ...
    Nw, nc_UKF_mc, EQM_final_UKF_mc, DtoP_final_UKF_mc, EQM_UKF_moy(end), DtoP_UKF_moy(end));

fclose(fid);

diagFile = fullfile(outDir, 'diagnostic_differences_EKF_UKF.txt');
fid = fopen(diagFile, 'w');

fprintf(fid, 'Diagnostic de differences EKF vs UKF\n');
fprintf(fid, '===================================\n\n');
fprintf(fid, 'Parametres : N = %d, Nw = %d, sigmau = %.6e, sigmav = %.6e\n', N, Nw, sigmau, sigmav);
fprintf(fid, 'UKF : alpha = %.6f, beta = %.6f, kappa = %.6f\n\n', alpha, beta, kappa);

fprintf(fid, 'Nw = 1\n');
fprintf(fid, 'max |EQM_UKF - EQM_EKF|       = %.10e\n', max(abs(EQM_UKF - EQM_EKF)));
fprintf(fid, 'max |DtoP_UKF - DtoP_EKF|     = %.10e\n', max(abs(DtoP_UKF - DtoP_EKF)));
fprintf(fid, 'max |traceP_UKF - traceP_EKF| = %.10e\n\n', max(abs(traceP_UKF - traceP_EKF)));

fprintf(fid, 'Monte-Carlo\n');
fprintf(fid, 'max |EQM_UKF_moy - EQM_EKF_moy|       = %.10e\n', max(abs(diff_EQM_MC)));
fprintf(fid, 'mean |EQM_UKF_moy - EQM_EKF_moy|      = %.10e\n', mean(abs(diff_EQM_MC)));
fprintf(fid, 'max |DtoP_UKF_moy - DtoP_EKF_moy|     = %.10e\n', max(abs(diff_DtoP_MC)));
fprintf(fid, 'mean |DtoP_UKF_moy - DtoP_EKF_moy|    = %.10e\n', mean(abs(diff_DtoP_MC)));
fprintf(fid, 'max |traceP_UKF_moy - traceP_EKF_moy| = %.10e\n', max(abs(diff_traceP_MC)));
fprintf(fid, 'mean |traceP_UKF_moy - traceP_EKF_moy|= %.10e\n', mean(abs(diff_traceP_MC)));

fclose(fid);


configFile = fullfile(outDir, 'parametres_valeurs_prof.txt');
fid = fopen(configFile, 'w');
fprintf(fid, 'Parametres alignes avec l''enonce de la prof\n');
fprintf(fid, '============================================\n');
fprintf(fid, 'Stations : Mk1=[-1500;1500], Mk2=[20000;20000], Mk3=[25000;500]\n');
fprintf(fid, 'Zinit = [1000;1000;0;0]\n');
fprintf(fid, 'Ts = %g, N = %d, S_seed = %d\n', Ts, N, S_seed);
fprintf(fid, 'sigmau = %g, sigmav = %.1e\n', sigmau, sigmav);
fprintf(fid, 'C1 : epsilon = %.4g, Nw = 1, z0 = [0;0;0;0]\n', epsilon);
fprintf(fid, 'C3 MC : epsilon = %.4g, Nw = %d\n', epsilon_MC, Nw);
fprintf(fid, 'UKF : alpha = %.4g, beta = %.4g, kappa = %.4g\n', alpha, beta, kappa);
fclose(fid);

fprintf('\nFigures, CSV et diagnostics sauvegardes dans :\n%s\n', outDir);
fprintf('Execution terminee.\n');

%% ========================================================================
%% Fonctions locales
%% ========================================================================

function [Zmat, Yf] = generate_traj_local(Zinit, sigmab, Ts, N, S_seed, Anchor_position)
% Genere la trajectoire vraie et les TOA non bruitees.
%
% Entrees :
%   Zinit           : etat initial [x; y; vx; vy]
%   sigmab          : variance du bruit d'acceleration
%   Ts              : periode d'echantillonnage
%   N               : nombre d'iterations
%   S_seed          : graine aleatoire
%   Anchor_position : positions des stations, taille 2 x Ns
%
% Sorties :
%   Zmat : trajectoire vraie, taille 4 x N
%   Yf   : TOA non bruitees, taille Ns x N

rng(S_seed);

A = [1 0 Ts 0;
     0 1 0 Ts;
     0 0 1 0;
     0 0 0 1];

B = [Ts^2/2 0;
     0 Ts^2/2;
     Ts 0;
     0 Ts];

c0 = 3e8;
Ns = size(Anchor_position, 2);

z = Zinit;
Zmat = zeros(4, N);
Yf = zeros(Ns, N);

for n = 1:N
    u = sqrt(sigmab) * randn(2, 1);
    z = A*z + B*u;

    Zmat(:, n) = z;

    for s = 1:Ns
        Yf(s, n) = norm(z(1:2) - Anchor_position(:, s)) / c0;
    end
end
end

function [Zmat_est, EQM, DtoP, normP_est] = EKF_track_traj_local( ...
    Zmat, Y_noisy, z_est, P_est, sigmau, sigmav, Ts, N, Anchor_position)
% Poursuite de trajectoire par EKF avec observations TOA.
%
% Difference principale avec UKF :
%   EKF linearise h(z) avec la Jacobienne H.
%
% Convention :
%   sigmau et sigmav sont des variances.

Ns = size(Anchor_position, 2);

A = [1 0 Ts 0;
     0 1 0 Ts;
     0 0 1 0;
     0 0 0 1];

B = [Ts^2/2 0;
     0 Ts^2/2;
     Ts 0;
     0 Ts];

Q = sigmau * (B * B');
R = sigmav * eye(Ns);
I4 = eye(4);

Zmat_est = zeros(4, N);
EQM = zeros(1, N);
DtoP = zeros(1, N);
normP_est = zeros(1, N);

err2_cum = 0;
pos_err_cum = 0;

for n = 1:N
    % 1) Prediction
    z_pred = A * z_est;
    P_pred = A * P_est * A' + Q;
    P_pred = make_symmetric_local(P_pred);

    % 2) Observation predite et Jacobienne
    y_pred = h_toa_local(z_pred, Anchor_position);
    H = jacobian_toa_local(z_pred, Anchor_position);

    % 3) Innovation et gain EKF
    S = H * P_pred * H' + R;
    S = make_symmetric_local(S) + 1e-18 * eye(Ns);

    K = P_pred * H' / S;
    innovation = Y_noisy(:, n) - y_pred;

    % 4) Correction
    z_est = z_pred + K * innovation;

    % Forme de Joseph pour ameliorer la stabilite numerique.
    P_est = (I4 - K*H) * P_pred * (I4 - K*H)' + K * R * K';
    P_est = make_symmetric_local(P_est);

    % 5) Stockage et performances
    Zmat_est(:, n) = z_est;

    err2_cum = err2_cum + norm(Zmat(:, n) - z_est)^2;
    pos_err_cum = pos_err_cum + norm(Zmat(1:2, n) - z_est(1:2));

    EQM(n) = err2_cum / n;
    DtoP(n) = pos_err_cum / n;
    normP_est(n) = trace(P_est);
end
end

function [Zmat_est, EQM, DtoP, normP_est] = UKF_track_traj_local( ...
    Zmat, Y_noisy, z_est, P_est, sigmau, sigmav, Ts, N, Anchor_position, alpha, beta, kappa)
% Poursuite de trajectoire par UKF avec observations TOA.
%
% Difference principale avec EKF :
%   UKF ne calcule pas la Jacobienne. Il propage des points sigma dans h(z).
%
% Convention :
%   sigmau et sigmav sont des variances.

L = length(z_est);
Ns = size(Anchor_position, 2);

A = [1 0 Ts 0;
     0 1 0 Ts;
     0 0 1 0;
     0 0 0 1];

B = [Ts^2/2 0;
     0 Ts^2/2;
     Ts 0;
     0 Ts];

Q = sigmau * (B * B');
R = sigmav * eye(Ns);

Zmat_est = zeros(4, N);
EQM = zeros(1, N);
DtoP = zeros(1, N);
normP_est = zeros(1, N);

err2_cum = 0;
pos_err_cum = 0;

for n = 1:N
    % 1) Points sigma autour de l'estimation courante
    [X, Wm, Wc] = sigma_points_local(z_est, P_est, alpha, beta, kappa);

    % 2) Prediction dynamique des points sigma
    X_pred = A * X;

    % 3) Moyenne predite
    z_pred = X_pred * Wm';

    % 4) Covariance predite
    P_pred = Q;
    for i = 1:(2*L + 1)
        dx = X_pred(:, i) - z_pred;
        P_pred = P_pred + Wc(i) * (dx * dx');
    end
    P_pred = make_symmetric_local(P_pred);

    % 5) Nouveaux points sigma autour de l'etat predit
    [X2, Wm, Wc] = sigma_points_local(z_pred, P_pred, alpha, beta, kappa);

    % 6) Projection non lineaire dans l'espace TOA
    Y_sigma = zeros(Ns, 2*L + 1);
    for i = 1:(2*L + 1)
        Y_sigma(:, i) = h_toa_local(X2(:, i), Anchor_position);
    end

    % 7) Mesure predite
    y_pred = Y_sigma * Wm';

    % 8) Covariance de l'innovation Pyy et covariance croisee Pxy
    Pyy = R;
    Pxy = zeros(L, Ns);

    for i = 1:(2*L + 1)
        dy = Y_sigma(:, i) - y_pred;
        dx = X2(:, i) - z_pred;

        Pyy = Pyy + Wc(i) * (dy * dy');
        Pxy = Pxy + Wc(i) * (dx * dy');
    end

    Pyy = make_symmetric_local(Pyy) + 1e-18 * eye(Ns);

    % 9) Gain UKF et correction
    K = Pxy / Pyy;
    innovation = Y_noisy(:, n) - y_pred;

    z_est = z_pred + K * innovation;
    P_est = P_pred - K * Pyy * K';
    P_est = make_symmetric_local(P_est);

    % 10) Stockage et performances
    Zmat_est(:, n) = z_est;

    err2_cum = err2_cum + norm(Zmat(:, n) - z_est)^2;
    pos_err_cum = pos_err_cum + norm(Zmat(1:2, n) - z_est(1:2));

    EQM(n) = err2_cum / n;
    DtoP(n) = pos_err_cum / n;
    normP_est(n) = trace(P_est);
end
end

function [X, Wm, Wc] = sigma_points_local(x, P, alpha, beta, kappa)
% Genere les points sigma et leurs poids.
%
% Pour un etat de dimension L, on genere 2L+1 points :
%   X0 = x
%   Xi = x + colonne_i(sqrt((L+lambda)P))
%   Xi = x - colonne_i(sqrt((L+lambda)P))
%
% Parametres :
%   alpha : dispersion des points sigma
%   beta  : information a priori sur la distribution, beta=2 pour gaussien
%   kappa : parametre secondaire

L = length(x);

lambda = alpha^2 * (L + kappa) - L;
c = L + lambda;

if c <= 0
    error('Parametres UKF invalides : L + lambda doit etre positif.');
end

Wm = zeros(1, 2*L + 1);
Wc = zeros(1, 2*L + 1);

Wm(1) = lambda / c;
Wc(1) = lambda / c + (1 - alpha^2 + beta);

Wm(2:end) = 1 / (2*c);
Wc(2:end) = 1 / (2*c);

P = make_symmetric_local(P);

% Stabilisation numerique de Cholesky
jitter = 1e-12;
success = false;

for tries = 1:10
    [S, flag] = chol(c * P + jitter * eye(L), 'lower');

    if flag == 0
        success = true;
        break;
    end

    jitter = jitter * 10;
end

if ~success
    % Dernier recours :
    % projection approximative vers une matrice semi-definie positive.
    [V, D] = eig(P);
    d = diag(D);
    d(d < 1e-12) = 1e-12;

    P = V * diag(d) * V';
    P = make_symmetric_local(P);

    S = chol(c * P + 1e-12 * eye(L), 'lower');
end

X = zeros(L, 2*L + 1);
X(:, 1) = x;

for i = 1:L
    X(:, i+1) = x + S(:, i);
    X(:, i+1+L) = x - S(:, i);
end
end

function y = h_toa_local(z, Anchor_position)
% Fonction d'observation TOA :
%   y_i = ||p - p_i|| / c0
%
% z contient :
%   z = [x; y; vx; vy]
%
% Anchor_position contient les stations en colonnes :
%   Anchor_position = [p1 p2 ... pNs]

c0 = 3e8;

Ns = size(Anchor_position, 2);
y = zeros(Ns, 1);

pos = z(1:2);

for s = 1:Ns
    y(s) = norm(pos - Anchor_position(:, s)) / c0;
end
end

function H = jacobian_toa_local(z, Anchor_position)
% Jacobienne de la fonction TOA h(z).
%
% Chaque ligne correspond a :
%   derivee de ||p - p_i|| / c0 par rapport a [x y vx vy]
%
% Les derivees par rapport aux vitesses sont nulles car la mesure TOA
% depend uniquement de la position instantanee.

c0 = 3e8;

Ns = size(Anchor_position, 2);
H = zeros(Ns, 4);

pos = z(1:2);

for s = 1:Ns
    diff = pos - Anchor_position(:, s);
    d = norm(diff);

    if d < 1e-9
        d = 1e-9;
    end

    H(s, 1:2) = (diff' / d) / c0;
    H(s, 3:4) = [0 0];
end
end

function M = make_symmetric_local(M)
% Force une matrice a etre symetrique numeriquement.
M = (M + M') / 2;
end

function [nc, EQM_final, DtoP_final] = final_perf_local(EQM, DtoP, eps_seuil)
% Detecte un indice de convergence nc puis calcule les performances finales.
%
% Critere utilise :
%   |EQM(k) - EQM(k-1)| < eps_seuil
%
% Si le critere n'est jamais atteint, nc reste egal a N.

N = length(EQM);
nc = N;

for k = 2:N
    if abs(EQM(k) - EQM(k-1)) < eps_seuil
        nc = k;
        break;
    end
end

EQM_final = mean(EQM(nc:N));
DtoP_final = mean(DtoP(nc:N));
end
