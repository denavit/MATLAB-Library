clear all; close all; clc;

fs = figureStyle('Display');

%% Define Cross Section
steel_shape = steel_shape_lookup('HSS12x12x1/4');

H  = steel_shape.Ht;
B  = steel_shape.B;
t  = steel_shape.tdes;
Fy = 46;
fc = 8;
units   = 'US';

section = RCFT(H,B,t,Fy,fc,units);
section.neglectLocalBuckling = true;

psd  = section.plasticStressDistributionObject;
aci  = section.strainCompatibilityAciObject;

%% Compute and Plot Two-Dimensional Interaction

[P_points,M_points] = section.sectionInteraction2d('strong','psd-acdbt','CompPos');

angle = 0;
numPoints = 50;
[P_cont,M_cont,~] = psd.interactionSweep(angle,numPoints);
[P_aci,M_aci,~]   = aci.interactionSweep(angle,numPoints);

hf = fs.figure(8,6);
ha = fs.axes();

plot(M_points/12,-P_points,'o--')
plot(M_cont/12,-P_cont)
plot(M_aci/12,-P_aci)
xlabel('Bending moment, M (kip-ft)')
zlabel('Axial compression, P (kips)')

xlim([0 320])
ylim([0 1500])

legend('PSD (discrete)','PSD (continuous)','ACI')

%% Compute and Plot Three-Dimensional Interaction 
numP        = 50;
numAngles   = 200;
numPoints_calc = 50;
numAngles_calc = 100;
results = aci.interaction3d_2(numP,numAngles,numPoints_calc,numAngles_calc);

hf = fs.figure(8,6);
ha = fs.axes();

surf(results.Mz(:,[1:end 1])/12,results.My(:,[1:end 1])/12,-results.P(:,[1:end 1]))
xlabel('Major-axis moment, Mz (kip-ft)')
ylabel('Minor-axis moment, My (kip-ft)')
zlabel('Axial compression, P (kips)')
view(150,26)

xlim([0 320])
ylim([0 320])
zlim([0 1500])
