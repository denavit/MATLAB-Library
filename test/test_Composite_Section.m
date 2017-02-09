clear all; close all; clc;

fs = figureStyle('Display');


%% Define Cross Section
steel_shape = steel_shape_lookup('HSS10.750x0.250');

D  = steel_shape.OD;
t  = steel_shape.tdes;
Fy = 42;
fc = 4;
units   = 'US';

section = CCFT(D,t,Fy,fc,units);
section.neglectLocalBuckling = true;


%% Plot Section
fs.picture(3,3)
section.plotSection;
axis_limits('Margin',0.1)


%% Compute Interaction
[P_points,M_points] = section.sectionInteraction2d('strong','psd-acdbt','CompPos');

psd  = section.plasticStressDistributionObject;
angle = 0;
numPoints = 50;
[P_cont,M_cont,~] = psd.interactionSweep(angle,numPoints);

aci  = section.strainCompatibilityAciObject;
angle = 0;
numPoints = 50;
[P_aci,M_aci,~] = aci.interactionSweep(angle,numPoints);


%% Plot Interaction 
fs.figure(5,5)
fs.axes;
plot(M_points,-P_points,'--ok','LineWidth',1)
plot(  M_cont,  -P_cont,  '-r','LineWidth',1)
plot(   M_aci,   -P_aci,  '-b','LineWidth',1)
xlabel('Bending Moment (kip-in.)')
ylabel('Axial Compression (kips)')

legend('PSD Points','PSD Cont.','ACI',...
    'Location','NE')

axis_limits('Margin',0.1)
axis_limits('Quadrant',1)
