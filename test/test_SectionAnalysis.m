clear all; close all; clc;

%% Read Section Definiton to Cell Array
def = strsplit(fileread('section.tcl'),'\n');
tag = 1;

%% Define SectionAnalysis Object
SA = SectionAnalysis(def,tag);
SA.echoOpenSeesOutput = true;

%% Discretization
SA.plotDiscretization();
SA.printMaterialInfo();

%% Material Stress Strain
matID = 1;
ss_results_tens = SA.getStressStrain([0 0.05],matID,'Steps',100);
ss_results_comp = SA.getStressStrain([0 -0.05],matID,'Steps',100);

figure
plot(ss_results_tens.disp,ss_results_tens.force)
hold all
plot(ss_results_comp.disp,ss_results_comp.force)

