clear all; close all; clc;

%% Material definition
% Units: kip, in, s
tag = 1;

Fy = 50;    % Yield strength
E0 = 29000; % Elastic modulus

bk = 0.01;  % Hardening ratio (kinematic hardening)
R0 = 20;    % Transition parameter from elastic to plastic branches
r1 = 0.90;  % Transition parameter from elastic to plastic branches
r2 = 0.15;  % Transition parameter from elastic to plastic branches

Fu = 65;    % Ultimate strength
Ru = 20;    % Transition to ultimate strength for kinematic hardening

matdef = sprintf('uniaxialMaterial Steel4 %d %g %g -kin %g %g %g %g -ult %g %g', ...
    tag, Fy, E0, bk, R0, r1, r2, Fu, Ru);

%% Analysis options
peakPoints  = [0 1 -1 2 -2 3 -3 4 -4 5 -5 0]*5*Fy/E0;
rateType    = 'StrainRate';
rateValue   = peakPoints(2)/100;

%% Run analysis
analysis = UniaxialMaterialAnalysis(matdef, tag);
analysis.deleteFilesAfterAnalysis = true;
analysis.echoOpenSeesOutput = false;

results_pos_env = analysis.runAnalysis([0 2*max(peakPoints)], rateType, rateValue);
results_neg_env = analysis.runAnalysis([0 2*min(peakPoints)], rateType, rateValue);
results         = analysis.runAnalysis(         peakPoints, rateType, rateValue);

%% Plot
figure
hold on
plot(results_pos_env.disp, results_pos_env.force, 'k--')
plot(results_neg_env.disp, results_neg_env.force, 'k--')
plot(results.disp, results.force)
xlabel('Strain (in/in)')
ylabel('Stress (ksi)')
