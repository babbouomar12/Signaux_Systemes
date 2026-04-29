clear; close all; clc;

%% ============================================================
%  COMPARAISON EKF vs UKF POUR L'ESTIMATION DU SOC D'UNE BATTERIE
%  + SIMULATION MONTE-CARLO POUR N_MC = 20, 50, 100
%
%  Objectif global :
%  -----------------
%  On veut estimer le SOC (State Of Charge) d'une batterie lithium-ion.
%  Le SOC représente le pourcentage réel de charge de la batterie.
%
%  Problème :
%  ----------
%  Le SOC n'est pas directement mesurable. En pratique, on mesure surtout :
%     - le courant I(k),
%     - la tension terminale V(k).
%
%  On utilise donc un filtre non linéaire pour reconstruire l'état caché :
%
%       x(k) = [ SOC(k) ; V_RC(k) ]
%
%  où :
%     - SOC(k)  : état de charge de la batterie,
%     - V_RC(k) : tension interne de polarisation, liée à la branche RC.
%
%  Deux filtres sont comparés :
%     - EKF : Extended Kalman Filter, basé sur la linéarisation par Jacobienne,
%     - UKF : Unscented Kalman Filter, basé sur les points sigma.
%
%  La simulation Monte-Carlo permet de vérifier que la comparaison EKF/UKF
%  ne dépend pas d'un seul tirage aléatoire du bruit.
%
%  Attention notation :
%  --------------------
%  Dans ce script :
%     - N = nombre d'échantillons temporels de la simulation.
%     - N_MC = nombre de réalisations Monte-Carlo.
%
%  Donc les valeurs 20, 50 et 100 concernent N_MC, pas N.
% =============================================================

%% ============================================================
%  0) Réglage général
% =============================================================

% On fixe la graine aléatoire pour que la simulation simple soit reproductible.
% Cela signifie que si on relance le script, on obtient les mêmes courbes
% pour la première simulation.
rng(7);

% Dossier où les figures seront sauvegardées automatiquement.
% Si le dossier n'existe pas, MATLAB le crée.
outDir = fullfile('figures', 'results');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

%% ============================================================
%  1) Paramètres temporels de la simulation
% =============================================================

% N est le nombre d'instants simulés.
% Ici, on simule 900 échantillons.
N  = 900;

% dt est la période d'échantillonnage en secondes.
% Ici, une mise à jour du modèle est effectuée chaque seconde.
dt = 1;

% t est le vecteur temps.
% Comme N = 900 et dt = 1, on obtient :
% t = [0, 1, 2, ..., 899].
t  = (0:N-1) * dt;

%% ============================================================
%  2) Paramètres physiques simplifiés de la batterie
% =============================================================

% Capacité nominale de la batterie en ampère-heure.
% 2.3 Ah correspond à une petite batterie lithium-ion.
Q_Ah = 2.3;

% Conversion de la capacité en coulombs.
% Pourquoi ? Parce que l'équation du SOC utilise une charge électrique.
% On sait que 1 Ah = 3600 C.
% Donc Q = 2.3 * 3600 = 8280 C.
Q = Q_Ah * 3600;

% Rendement coulombique.
% eta = 1 signifie qu'on suppose un rendement parfait.
% Toute la charge consommée est donc retirée de la batterie.
eta = 1.0;

% Résistance ohmique interne de la batterie.
% Elle provoque une chute de tension instantanée R0*I.
R0 = 0.045;

% Résistance de polarisation de la branche RC.
% Elle modélise une partie du comportement dynamique interne de la batterie.
Rp = 0.025;

% Capacité de polarisation de la branche RC.
% Avec Rp, elle fixe la vitesse de réponse de V_RC.
Cp = 2500;

% Coefficient discret de la branche RC.
% Le modèle continu d'une branche RC donne une réponse exponentielle.
% Après discrétisation :
%
%       V_RC(k+1) = a V_RC(k) + Rp(1-a) I(k)
%
% avec :
%
%       a = exp(-dt/(Rp*Cp))
%
a = exp(-dt/(Rp*Cp));

%% ============================================================
%  3) Bruits de modèle et de mesure
% =============================================================

% sigmaV est l'écart-type du bruit de mesure de tension.
% Ici sigmaV = 0.025 V, donc environ 25 mV.
% La mesure de tension sera :
%
%       V_meas(k) = V_true(k) + sigmaV * randn
%
sigmaV = 0.025;

% Q_filter est la covariance du bruit de modèle utilisée par EKF et UKF.
% Elle représente ce que les filtres CROIENT sur l'imperfection du modèle.
%
% L'état est :
%
%       x = [ SOC ; V_RC ]
%
% Donc Q_filter est une matrice 2x2.
%
% diag([2e-7, 5e-6]) donne :
%
%       [ 2e-7    0    ]
%       [  0     5e-6  ]
%
% Le premier terme concerne l'incertitude sur SOC.
% Le deuxième terme concerne l'incertitude sur V_RC.
Q_filter = diag([2e-7, 5e-6]);

% R_filter est la variance du bruit de mesure utilisée par les filtres.
% Le filtre de Kalman utilise une variance, pas directement l'écart-type.
%
%       R = sigmaV^2
%
R_filter = sigmaV^2;

% Q_true est la covariance du bruit réellement injecté dans la simulation
% de la vraie batterie.
%
% On distingue Q_true et Q_filter car, dans la réalité, le filtre ne connaît
% jamais exactement le vrai bruit du système.
Q_true = diag([5e-8, 2e-6]);

%% ============================================================
%  4) Courant de décharge I(k)
% =============================================================

% Convention :
% ------------
% I(k) > 0 signifie que la batterie se décharge.
%
% Cela est cohérent avec l'équation :
%
%       SOC(k+1) = SOC(k) - eta*dt*I(k)/Q
%
% Si I(k) est positif, alors le SOC diminue.
%
% On crée un profil de courant variable pour rendre la simulation plus
% réaliste et plus intéressante qu'une simple décharge constante.
I = create_current_profile(t, N);

%% ============================================================
%  5) Initialisation des filtres
% =============================================================

% Etat initial utilisé par EKF et UKF.
% On choisit volontairement une mauvaise estimation initiale pour tester
% la capacité des filtres à converger.
%
% La vraie batterie commencera à SOC = 95 %.
% Mais les filtres croient au départ que SOC = 60 %.
x0_est = [0.60; 0.08];

% Covariance initiale.
% Elle indique l'incertitude initiale du filtre.
%
% 0.20^2 signifie que l'on accepte une forte incertitude initiale sur le SOC.
% 0.08^2 signifie une incertitude plus faible sur V_RC.
P0 = diag([0.20^2, 0.08^2]);

%% ============================================================
%  6) Paramètres UKF
% =============================================================

% Dimension de l'état :
% x = [SOC ; V_RC], donc n = 2.
n = 2;

% alpha règle l'écartement des points sigma autour de la moyenne.
% Plus alpha est grand, plus les sigma-points s'éloignent de l'estimation.
alpha = 0.7;

% beta incorpore une information sur la distribution.
% beta = 2 est le choix classique pour une distribution gaussienne.
beta = 2;

% kappa est un paramètre secondaire de mise à l'échelle.
% kappa = 0 est un choix courant en simulation.
kappa = 0;

%% ============================================================
%  7) Regroupement des paramètres dans une structure
% =============================================================

% On range les paramètres dans une structure "params" pour les passer
% facilement aux fonctions locales, notamment à run_one_soc_trial.
params.N = N;
params.dt = dt;

params.Q = Q;
params.eta = eta;
params.R0 = R0;
params.Rp = Rp;
params.Cp = Cp;
params.a = a;

params.sigmaV = sigmaV;
params.Q_filter = Q_filter;
params.R_filter = R_filter;
params.Q_true = Q_true;

params.x0_est = x0_est;
params.P0 = P0;

params.alpha = alpha;
params.beta = beta;
params.kappa = kappa;

%% ============================================================
%  8) Simulation simple : une seule réalisation du bruit
% =============================================================

% Cette simulation permet de visualiser une trajectoire particulière :
% - SOC réel,
% - SOC estimé EKF,
% - SOC estimé UKF.
%
% Le troisième argument, ici 7, est la graine aléatoire utilisée pour générer
% les bruits de cette réalisation.
data = run_one_soc_trial(params, I, 7);

% On récupère les grandeurs utiles dans des variables plus lisibles.
soc_true = data.soc_true;
soc_ekf  = data.soc_ekf;
soc_ukf  = data.soc_ukf;

% Erreur d'estimation :
% erreur = SOC réel - SOC estimé.
err_ekf = soc_true - soc_ekf;
err_ukf = soc_true - soc_ukf;

% Affichage des RMSE dans la console.
fprintf('\n================ Résultats simulation simple ================\n');
fprintf('RMSE EKF sur le SOC : %.5f soit %.3f %%\n', data.rmse_ekf, data.rmse_ekf*100);
fprintf('RMSE UKF sur le SOC : %.5f soit %.3f %%\n', data.rmse_ukf, data.rmse_ukf*100);
fprintf('Gain UKF vs EKF     : %.2f %%\n', 100*(data.rmse_ekf-data.rmse_ukf)/data.rmse_ekf);
fprintf('=============================================================\n\n');

%% ============================================================
%  9) Figures de la simulation simple
% =============================================================

% Figure du courant de décharge.
figure('Name','Courant de décharge');
plot(t, I, 'LineWidth', 1.4);
grid on;
xlabel('Temps [s]');
ylabel('Courant I [A]');
title('Profil de courant de décharge');
save_current_figure(outDir, 'fig_09_courant_decharge.png');

% Figure principale : SOC réel vs SOC estimé par EKF et UKF.
figure('Name','SOC réel vs estimations EKF et UKF');
plot(t, soc_true*100, 'k', 'LineWidth', 2); hold on;
plot(t, soc_ekf*100, '--', 'LineWidth', 1.5);
plot(t, soc_ukf*100, '-.', 'LineWidth', 1.5);
grid on;
xlabel('Temps [s]');
ylabel('SOC [%]');
title('Estimation du SOC : EKF vs UKF');
legend('SOC réel', 'SOC estimé EKF', 'SOC estimé UKF', 'Location', 'best');
save_current_figure(outDir, 'fig_10_SOC_EKF_vs_UKF.png');

% Figure des erreurs d'estimation.
figure('Name','Erreur SOC');
plot(t, err_ekf*100, '--', 'LineWidth', 1.5); hold on;
plot(t, err_ukf*100, '-.', 'LineWidth', 1.5);
plot([t(1), t(end)], [0, 0], 'k', 'LineWidth', 1.0);
grid on;
xlabel('Temps [s]');
ylabel('Erreur SOC [%]');
title('Erreur d''estimation du SOC');
legend('Erreur EKF', 'Erreur UKF', 'Erreur nulle', 'Location', 'best');
save_current_figure(outDir, 'fig_11_erreur_SOC_EKF_vs_UKF.png');

% Figure de tension.
% On compare :
% - la tension vraie,
% - la tension mesurée bruitée,
% - la tension reconstruite à partir de l'état estimé par EKF,
% - la tension reconstruite à partir de l'état estimé par UKF.
figure('Name','Tension mesurée et tension reconstruite');
plot(t, data.V_meas, '.', 'MarkerSize', 5); hold on;
plot(t, data.V_true, 'k', 'LineWidth', 2);
plot(t, data.V_ekf, '--', 'LineWidth', 1.4);
plot(t, data.V_ukf, '-.', 'LineWidth', 1.4);
grid on;
xlabel('Temps [s]');
ylabel('Tension [V]');
title('Tension mesurée et tension estimée');
legend('Tension mesurée bruitée', 'Tension réelle', ...
       'Tension reconstruite EKF', 'Tension reconstruite UKF', ...
       'Location', 'best');
save_current_figure(outDir, 'fig_12_tension_mesuree_reconstruite.png');

% Figure RMSE en barres.
figure('Name','Comparaison RMSE');
bar([data.rmse_ekf, data.rmse_ukf]*100);
grid on;
set(gca, 'XTickLabel', {'EKF', 'UKF'});
ylabel('RMSE SOC [%]');
title('Comparaison quantitative EKF vs UKF');
save_current_figure(outDir, 'fig_13_RMSE_SOC_EKF_vs_UKF.png');

% Figure de la courbe OCV(SOC).
% Cette figure justifie la non-linéarité du modèle de mesure.
figure('Name','Courbe OCV(SOC)');
soc_grid = linspace(0,1,500);
plot(soc_grid*100, OCV(soc_grid), 'LineWidth', 2);
grid on;
xlabel('SOC [%]');
ylabel('OCV [V]');
title('Courbe non linéaire OCV(SOC) utilisée dans la simulation');
save_current_figure(outDir, 'fig_14_courbe_OCV_SOC.png');

%% ============================================================
%  10) Simulation Monte-Carlo pour N_MC = 20, 50 et 100
% =============================================================

% Ici, on veut vérifier la robustesse statistique de la comparaison EKF/UKF.
%
% Pour chaque valeur de N_MC :
%   - on répète la simulation N_MC fois,
%   - à chaque fois, on régénère les bruits,
%   - EKF et UKF utilisent les mêmes mesures bruitées dans une même réalisation,
%   - on calcule la RMSE pour chaque réalisation,
%   - on moyenne les RMSE obtenues.
%
% On teste trois nombres de réalisations Monte-Carlo.
Nmc_list = [20, 50, 100];

% Préallocation des résultats.
% Chaque tableau contient une valeur pour N_MC=20, une pour 50, une pour 100.
mean_rmse_ekf_list = zeros(size(Nmc_list));
mean_rmse_ukf_list = zeros(size(Nmc_list));

std_rmse_ekf_list = zeros(size(Nmc_list));
std_rmse_ukf_list = zeros(size(Nmc_list));

% Pour afficher aussi les résultats détaillés de chaque réalisation, on garde
% les RMSE dans des cellules. Chaque cellule correspond à une valeur de N_MC.
rmse_ekf_all = cell(size(Nmc_list));
rmse_ukf_all = cell(size(Nmc_list));

% Boucle externe : on parcourt les valeurs 20, 50 et 100.
for idx = 1:length(Nmc_list)

    % Nombre de réalisations Monte-Carlo courant.
    Nmc = Nmc_list(idx);

    % Vecteurs qui stockent les RMSE pour chaque réalisation.
    rmse_ekf_mc = zeros(1, Nmc);
    rmse_ukf_mc = zeros(1, Nmc);

    % Boucle Monte-Carlo.
    for mc = 1:Nmc

        % Graine différente pour chaque réalisation.
        % 1000*idx permet de séparer les séries de tirages pour N_MC=20,50,100.
        seed_mc = 1000*idx + mc;

        % On lance une simulation complète : vraie batterie + mesures bruitées + EKF + UKF.
        data_mc = run_one_soc_trial(params, I, seed_mc);

        % On stocke les RMSE de cette réalisation.
        rmse_ekf_mc(mc) = data_mc.rmse_ekf;
        rmse_ukf_mc(mc) = data_mc.rmse_ukf;
    end

    % On garde les RMSE brutes dans des cellules pour d'éventuelles figures.
    rmse_ekf_all{idx} = rmse_ekf_mc;
    rmse_ukf_all{idx} = rmse_ukf_mc;

    % Moyenne des RMSE sur les Nmc réalisations.
    mean_rmse_ekf_list(idx) = mean(rmse_ekf_mc);
    mean_rmse_ukf_list(idx) = mean(rmse_ukf_mc);

    % Ecart-type des RMSE.
    % Il donne une idée de la dispersion des performances.
    std_rmse_ekf_list(idx) = std(rmse_ekf_mc);
    std_rmse_ukf_list(idx) = std(rmse_ukf_mc);
end

%% Affichage console des résultats Monte-Carlo

fprintf('\n================ Monte-Carlo multi N_MC ================\n');
for idx = 1:length(Nmc_list)
    fprintf('N_MC = %d réalisations\n', Nmc_list(idx));
    fprintf('  RMSE moyenne EKF = %.5f soit %.3f %%\n', ...
        mean_rmse_ekf_list(idx), mean_rmse_ekf_list(idx)*100);
    fprintf('  RMSE moyenne UKF = %.5f soit %.3f %%\n', ...
        mean_rmse_ukf_list(idx), mean_rmse_ukf_list(idx)*100);
    fprintf('  Ecart-type EKF   = %.5f soit %.3f %%\n', ...
        std_rmse_ekf_list(idx), std_rmse_ekf_list(idx)*100);
    fprintf('  Ecart-type UKF   = %.5f soit %.3f %%\n', ...
        std_rmse_ukf_list(idx), std_rmse_ukf_list(idx)*100);
    fprintf('  Gain moyen UKF   = %.2f %%\n\n', ...
        100*(mean_rmse_ekf_list(idx)-mean_rmse_ukf_list(idx))/mean_rmse_ekf_list(idx));
end
fprintf('========================================================\n\n');

%% ============================================================
%  11) Figures Monte-Carlo pour N_MC = 20, 50, 100
% =============================================================

% Figure 1 :
% RMSE moyenne en fonction du nombre de réalisations Monte-Carlo.
figure('Name','Monte-Carlo : RMSE moyenne selon N_MC');
plot(Nmc_list, mean_rmse_ekf_list*100, 'o-', 'LineWidth', 1.7); hold on;
plot(Nmc_list, mean_rmse_ukf_list*100, 's-', 'LineWidth', 1.7);
grid on;
xlabel('Nombre de réalisations Monte-Carlo N_{MC}');
ylabel('RMSE moyenne du SOC [%]');
title('Influence de N_{MC} sur la RMSE moyenne');
legend('EKF', 'UKF', 'Location', 'best');
save_current_figure(outDir, 'fig_18_MC_multi_Nmc_RMSE_moyenne.png');

% Figure 2 :
% RMSE moyenne avec barres d'erreur.
% Les barres d'erreur correspondent à l'écart-type des RMSE.
figure('Name','Monte-Carlo : RMSE moyenne avec écart-type');
errorbar(Nmc_list, mean_rmse_ekf_list*100, std_rmse_ekf_list*100, ...
         'o-', 'LineWidth', 1.7); hold on;
errorbar(Nmc_list, mean_rmse_ukf_list*100, std_rmse_ukf_list*100, ...
         's-', 'LineWidth', 1.7);
grid on;
xlabel('Nombre de réalisations Monte-Carlo N_{MC}');
ylabel('RMSE du SOC [%]');
title('RMSE moyenne et dispersion selon N_{MC}');
legend('EKF', 'UKF', 'Location', 'best');
save_current_figure(outDir, 'fig_19_MC_multi_Nmc_RMSE_ecart_type.png');

% Figure 3 :
% Affichage des RMSE de chaque réalisation pour N_MC = 100.
% Cela montre la dispersion réalisation par réalisation.
idx100 = find(Nmc_list == 100, 1);
if ~isempty(idx100)
    figure('Name','Monte-Carlo : RMSE par réalisation pour N_MC = 100');
    plot(1:100, rmse_ekf_all{idx100}*100, 'o-', 'LineWidth', 1.2); hold on;
    plot(1:100, rmse_ukf_all{idx100}*100, 's-', 'LineWidth', 1.2);
    grid on;
    xlabel('Réalisation Monte-Carlo');
    ylabel('RMSE SOC [%]');
    title('RMSE du SOC pour chaque réalisation, N_{MC}=100');
    legend('EKF', 'UKF', 'Location', 'best');
    save_current_figure(outDir, 'fig_20_MC_RMSE_par_realisation_Nmc100.png');
end

disp('Toutes les figures ont été sauvegardées dans le dossier figures/results/.');

%% ============================================================
%  FONCTIONS LOCALES
% =============================================================

function I = create_current_profile(t, N)
    %% Création du courant de décharge
    %
    % Cette fonction crée un courant artificiel mais utile pédagogiquement.
    % On veut éviter un courant constant, car un courant constant produit une
    % dynamique trop simple.
    %
    % Le courant contient :
    %   - une composante moyenne,
    %   - une variation lente,
    %   - une variation rapide,
    %   - des pulses périodiques,
    %   - un pulse plus fort entre k=520 et k=660.
    %
    % Convention :
    %   I > 0 signifie décharge.

    % Courant de base :
    %
    % 0.7 A : valeur moyenne de décharge.
    % 0.35*sin(...) : variation lente de période 130 s.
    % 0.20*sin(...) : variation plus rapide de période 37 s.
    I = 0.7 ...
        + 0.35*sin(2*pi*t/130) ...
        + 0.20*sin(2*pi*t/37);

    % Ajout de pulses.
    for k = 1:N

        % floor(k/80) divise le temps en blocs de 80 échantillons.
        % mod(...,2) alterne entre 0 et 1.
        % Lorsque le résultat vaut 1, on ajoute 0.9 A.
        if mod(floor(k/80), 2) == 1
            I(k) = I(k) + 0.9;
        end

        % Gros pulse entre k=520 et k=660.
        % Il simule une phase de forte sollicitation de la batterie.
        if k > 520 && k < 660
            I(k) = I(k) + 1.2;
        end
    end

    % On impose un courant minimum positif.
    % Cela évite de passer en régime de charge.
    I = max(I, 0.05);
end

function data = run_one_soc_trial(params, I, seed)
    %% Une réalisation complète du problème batterie
    %
    % Cette fonction effectue une simulation complète :
    %
    %   1) simulation du système réel,
    %   2) génération des mesures bruitées,
    %   3) estimation par EKF,
    %   4) estimation par UKF,
    %   5) calcul des erreurs et de la RMSE.
    %
    % Le paramètre seed fixe la graine aléatoire de cette réalisation.

    rng(seed);

    %% Récupération des paramètres
    %
    % On extrait les champs de la structure params pour rendre le code
    % plus lisible dans la suite.

    N = params.N;
    dt = params.dt;

    Q = params.Q;
    eta = params.eta;
    a = params.a;
    Rp = params.Rp;
    R0 = params.R0;

    sigmaV = params.sigmaV;
    Q_true = params.Q_true;
    Q_filter = params.Q_filter;
    R_filter = params.R_filter;

    x0_est = params.x0_est;
    P0 = params.P0;

    alpha = params.alpha;
    beta  = params.beta;
    kappa = params.kappa;

    %% ========================================================
    %  A) Simulation du système réel
    % =========================================================

    % x_true stocke l'état réel de la batterie.
    %
    % Ligne 1 : SOC réel.
    % Ligne 2 : V_RC réel.
    x_true = zeros(2, N);

    % Etat réel initial :
    % SOC = 95 %, V_RC = 0.
    x_true(:,1) = [0.95; 0.0];

    % V_true : tension idéale sans bruit.
    % V_meas : tension mesurée avec bruit.
    V_true = zeros(1, N);
    V_meas = zeros(1, N);

    % Propagation du vrai système.
    for k = 2:N

        % Application du modèle dynamique de la batterie.
        % On calcule x_true(k) à partir de x_true(k-1) et du courant I(k-1).
        x_true(:,k) = f_battery(x_true(:,k-1), I(k-1), dt, Q, eta, a, Rp);

        % Ajout d'un bruit de modèle.
        %
        % randn(2,1) génère deux bruits gaussiens standards.
        % chol(Q_true,'lower') permet de leur donner la covariance Q_true.
        x_true(:,k) = x_true(:,k) + chol(Q_true, 'lower') * randn(2,1);

        % Saturation physique du SOC.
        % Le SOC doit rester dans un intervalle réaliste.
        x_true(1,k) = clip(x_true(1,k), 0.01, 0.99);
    end

    % Génération des tensions.
    for k = 1:N

        % Tension vraie, sans bruit de capteur.
        V_true(k) = h_battery(x_true(:,k), I(k), R0);

        % Tension mesurée, avec bruit de mesure.
        % sigmaV est l'écart-type du bruit.
        V_meas(k) = V_true(k) + sigmaV * randn;
    end

    %% ========================================================
    %  B) EKF : Extended Kalman Filter
    % =========================================================

    % x_ekf stocke les estimations EKF de l'état.
    x_ekf = zeros(2, N);

    % Initialisation EKF.
    x_ekf(:,1) = x0_est;
    P_ekf = P0;

    for k = 2:N

        % ---------------- Prédiction de l'état ----------------
        %
        % On prédit l'état suivant avec le modèle dynamique.
        x_pred = f_battery(x_ekf(:,k-1), I(k-1), dt, Q, eta, a, Rp);

        % Jacobienne du modèle dynamique.
        %
        % SOC(k+1) dépend linéairement de SOC(k) avec dérivée 1.
        % V_RC(k+1) dépend de V_RC(k) avec dérivée a.
        F = [1, 0;
             0, a];

        % Prédiction de la covariance.
        %
        % P_pred représente l'incertitude après prédiction.
        P_pred = F * P_ekf * F' + Q_filter;

        % ---------------- Correction par la mesure -------------
        %
        % Modèle de mesure :
        %
        %    V = OCV(SOC) - R0*I - V_RC
        %
        % Jacobienne H = [dOCV/dSOC, -1].
        H = [dOCV_dSOC(x_pred(1)), -1];

        % Tension prédite par l'état prédit.
        z_pred = h_battery(x_pred, I(k), R0);

        % Covariance de l'innovation.
        %
        % Elle représente l'incertitude sur l'écart mesure-prédiction.
        S = H * P_pred * H' + R_filter;

        % Gain de Kalman.
        %
        % Il détermine combien la mesure doit corriger la prédiction.
        K = P_pred * H' / S;

        % Innovation :
        %
        % différence entre la tension mesurée et la tension prédite.
        innovation = V_meas(k) - z_pred;

        % Correction de l'état.
        x_upd = x_pred + K * innovation;

        % Correction de la covariance avec la forme de Joseph.
        % Cette forme est plus stable numériquement.
        eye2 = eye(2);
        P_upd = (eye2 - K*H) * P_pred * (eye2 - K*H)' + K * R_filter * K';

        % Saturation physique du SOC estimé.
        x_upd(1) = clip(x_upd(1), 0.01, 0.99);

        % Stockage.
        x_ekf(:,k) = x_upd;

        % On force la covariance à rester symétrique.
        P_ekf = (P_upd + P_upd') / 2;
    end

    %% ========================================================
    %  C) UKF : Unscented Kalman Filter
    % =========================================================

    % x_ukf stocke les estimations UKF.
    x_ukf = zeros(2, N);

    % Initialisation UKF identique à l'EKF pour une comparaison équitable.
    x_ukf(:,1) = x0_est;
    P_ukf = P0;

    % Dimension de l'état.
    n = 2;

    % Paramètre lambda de la transformation Unscented.
    lambda = alpha^2 * (n + kappa) - n;

    % Poids pour la moyenne et la covariance.
    Wm = zeros(1, 2*n + 1);
    Wc = zeros(1, 2*n + 1);

    % Poids du point central.
    Wm(1) = lambda / (n + lambda);
    Wc(1) = lambda / (n + lambda) + (1 - alpha^2 + beta);

    % Poids des autres points sigma.
    for i = 2:(2*n + 1)
        Wm(i) = 1 / (2*(n + lambda));
        Wc(i) = 1 / (2*(n + lambda));
    end

    for k = 2:N

        % ---------------- Génération des sigma-points ----------
        %
        % Les sigma-points représentent l'incertitude autour de x_ukf(:,k-1).
        Xsig = sigma_points(x_ukf(:,k-1), P_ukf, lambda);

        % ---------------- Prédiction des sigma-points ----------
        %
        % Chaque sigma-point est propagé dans le modèle dynamique.
        Xpred = zeros(n, 2*n + 1);

        for i = 1:(2*n + 1)
            Xpred(:,i) = f_battery(Xsig(:,i), I(k-1), dt, Q, eta, a, Rp);
        end

        % Moyenne prédite.
        x_pred = Xpred * Wm';

        % Covariance prédite.
        P_pred = Q_filter;

        for i = 1:(2*n + 1)
            dx = Xpred(:,i) - x_pred;
            P_pred = P_pred + Wc(i) * (dx * dx');
        end

        % Symétrisation numérique.
        P_pred = (P_pred + P_pred') / 2;

        % ---------------- Projection dans l'espace de mesure ---
        %
        % Chaque sigma-point prédit est transformé en tension prédite.
        Zsig = zeros(1, 2*n + 1);

        for i = 1:(2*n + 1)
            Zsig(i) = h_battery(Xpred(:,i), I(k), R0);
        end

        % Tension moyenne prédite.
        z_pred = Zsig * Wm';

        % Variance de l'innovation.
        S = R_filter;

        % Covariance croisée état-mesure.
        % Elle permet de convertir une erreur de tension en correction d'état.
        Pxz = zeros(n,1);

        for i = 1:(2*n + 1)
            dz = Zsig(i) - z_pred;
            dx = Xpred(:,i) - x_pred;

            S   = S   + Wc(i) * dz^2;
            Pxz = Pxz + Wc(i) * dx * dz;
        end

        % ---------------- Correction UKF ----------------------
        %
        % Gain de Kalman UKF.
        K = Pxz / S;

        % Innovation mesure - prédiction.
        innovation = V_meas(k) - z_pred;

        % Correction de l'état.
        x_upd = x_pred + K * innovation;

        % Correction de la covariance.
        P_upd = P_pred - K * S * K';

        % Saturation physique du SOC.
        x_upd(1) = clip(x_upd(1), 0.01, 0.99);

        % Stockage.
        x_ukf(:,k) = x_upd;

        % Symétrisation.
        P_ukf = (P_upd + P_upd') / 2;
    end

    %% ========================================================
    %  D) Reconstruction de la tension estimée
    % =========================================================

    % Ces tensions sont reconstruites à partir des états estimés.
    % Elles servent seulement à vérifier si l'état estimé explique bien
    % la tension mesurée.
    V_ekf = zeros(1, N);
    V_ukf = zeros(1, N);

    for k = 1:N
        V_ekf(k) = h_battery(x_ekf(:,k), I(k), R0);
        V_ukf(k) = h_battery(x_ukf(:,k), I(k), R0);
    end

    %% ========================================================
    %  E) Calcul des erreurs et RMSE
    % =========================================================

    % Extraction des SOC.
    soc_true = x_true(1,:);
    soc_ekf  = x_ekf(1,:);
    soc_ukf  = x_ukf(1,:);

    % Erreurs.
    err_ekf = soc_true - soc_ekf;
    err_ukf = soc_true - soc_ukf;

    % RMSE.
    % Comme SOC est entre 0 et 1, on multiplie souvent par 100 pour l'exprimer en %.
    rmse_ekf = sqrt(mean(err_ekf.^2));
    rmse_ukf = sqrt(mean(err_ukf.^2));

    %% ========================================================
    %  F) Sortie
    % =========================================================

    data.x_true = x_true;
    data.x_ekf = x_ekf;
    data.x_ukf = x_ukf;

    data.soc_true = soc_true;
    data.soc_ekf = soc_ekf;
    data.soc_ukf = soc_ukf;

    data.V_true = V_true;
    data.V_meas = V_meas;
    data.V_ekf = V_ekf;
    data.V_ukf = V_ukf;

    data.rmse_ekf = rmse_ekf;
    data.rmse_ukf = rmse_ukf;
end

function x_next = f_battery(x, I, dt, Q, eta, a, Rp)
    %% Modèle dynamique de la batterie
    %
    % Entrée :
    %   x = [SOC ; V_RC]
    %
    % Equations :
    %
    %   SOC(k+1) = SOC(k) - eta*dt*I(k)/Q
    %
    %   V_RC(k+1) = a*V_RC(k) + Rp*(1-a)*I(k)
    %
    % Si I est positif, la batterie se décharge et le SOC diminue.

    % Extraction du SOC.
    soc = x(1);

    % Extraction de la tension de polarisation.
    vrc = x(2);

    % Evolution du SOC.
    soc_next = soc - eta * dt * I / Q;

    % Evolution de la tension de polarisation.
    vrc_next = a * vrc + Rp * (1 - a) * I;

    % Saturation du SOC dans un intervalle physique.
    soc_next = clip(soc_next, 0.01, 0.99);

    % Reconstruction du vecteur d'état suivant.
    x_next = [soc_next; vrc_next];
end

function V = h_battery(x, I, R0)
    %% Modèle de mesure de la tension
    %
    % La tension terminale est modélisée par :
    %
    %   V = OCV(SOC) - R0*I - V_RC
    %
    % où :
    %   OCV(SOC) est la tension à vide,
    %   R0*I est la chute ohmique,
    %   V_RC est la chute dynamique interne.

    % Extraction du SOC.
    soc = x(1);

    % Extraction de V_RC.
    vrc = x(2);

    % Calcul de la tension terminale.
    V = OCV(soc) - R0 * I - vrc;
end

function y = OCV(soc)
    %% Courbe OCV(SOC)
    %
    % Cette fonction donne la tension à vide de la batterie en fonction du SOC.
    % Dans une vraie batterie, cette courbe est obtenue expérimentalement.
    %
    % Ici, on utilise une courbe artificielle mais non linéaire pour montrer
    % l'intérêt de comparer EKF et UKF.

    % On évite les valeurs exactement 0 ou 1 pour limiter les problèmes numériques.
    soc = clip(soc, 0.001, 0.999);

    % Courbe non linéaire :
    % - 3.15 : tension de base,
    % - 0.55*soc : composante globalement croissante,
    % - tanh autour de 0.15 : transition rapide à bas SOC,
    % - tanh autour de 0.85 : transition rapide à haut SOC,
    % - sinus : petite ondulation pour renforcer la non-linéarité.
    y = 3.15 ...
        + 0.55*soc ...
        + 0.18*tanh(20*(soc - 0.15)) ...
        + 0.15*tanh(18*(soc - 0.85)) ...
        + 0.015*sin(8*pi*soc);
end

function dy = dOCV_dSOC(soc)
    %% Dérivée analytique de OCV(SOC)
    %
    % Cette dérivée est utilisée uniquement par l'EKF.
    % L'UKF n'en a pas besoin, car il utilise les points sigma.

    % Saturation numérique.
    soc = clip(soc, 0.001, 0.999);

    % Dérivée de 0.55*soc.
    term1 = 0.55;

    % Arguments internes des fonctions tanh.
    z1 = 20*(soc - 0.15);
    z2 = 18*(soc - 0.85);

    % Dérivée de tanh(u) :
    %
    %   d/dx tanh(u(x)) = u'(x) / cosh(u(x))^2
    %
    term2 = 0.18 * 20 * (1 ./ cosh(z1)).^2;
    term3 = 0.15 * 18 * (1 ./ cosh(z2)).^2;

    % Dérivée du sinus.
    term4 = 0.015 * 8*pi * cos(8*pi*soc);

    % Dérivée totale.
    dy = term1 + term2 + term3 + term4;
end

function Xsig = sigma_points(x, P, lambda)
    %% Génération des sigma-points de l'UKF
    %
    % Pour un état de dimension n, l'UKF crée 2n+1 sigma-points :
    %
    %   X0       = x
    %   Xi       = x + colonne_i(sqrt((n+lambda)P))
    %   X{i+n}   = x - colonne_i(sqrt((n+lambda)P))
    %
    % Ces points représentent la moyenne et la covariance de l'état.

    % Dimension de l'état.
    n = numel(x);

    % On force P à être symétrique.
    % Cela évite des problèmes numériques dans chol.
    P = (P + P') / 2;

    % Petit terme ajouté à la diagonale si nécessaire.
    % Il aide lorsque P est presque définie positive mais pas parfaitement.
    jitter = 1e-10;
    success = false;

    % On tente la décomposition de Cholesky.
    for attempt = 1:8

        % A est la racine matricielle de (n+lambda)P.
        [A, p] = chol((n + lambda) * (P + jitter*eye(n)), 'lower');

        % Si p == 0, la décomposition a réussi.
        if p == 0
            success = true;
            break;
        end

        % Sinon, on augmente légèrement le jitter.
        jitter = jitter * 10;
    end

    % Si la décomposition échoue encore, on arrête le programme.
    if ~success
        error('La matrice de covariance n''est pas définie positive.');
    end

    % Préallocation de la matrice des sigma-points.
    Xsig = zeros(n, 2*n + 1);

    % Premier point sigma : la moyenne elle-même.
    Xsig(:,1) = x;

    % Points positifs et négatifs autour de la moyenne.
    for i = 1:n
        Xsig(:,i+1)   = x + A(:,i);
        Xsig(:,i+1+n) = x - A(:,i);
    end
end

function y = clip(x, xmin, xmax)
    %% Saturation d'une valeur
    %
    % Cette fonction force x à rester entre xmin et xmax.
    %
    % Si x < xmin, elle renvoie xmin.
    % Si x > xmax, elle renvoie xmax.
    % Sinon, elle renvoie x.

    y = min(max(x, xmin), xmax);
end

function save_current_figure(outDir, fileName)
    %% Sauvegarde robuste d'une figure
    %
    % exportgraphics donne une bonne qualité d'image sur les versions récentes.
    % saveas est utilisé comme solution de secours.

    fullName = fullfile(outDir, fileName);

    try
        exportgraphics(gcf, fullName, 'Resolution', 300);
    catch
        saveas(gcf, fullName);
    end
end
