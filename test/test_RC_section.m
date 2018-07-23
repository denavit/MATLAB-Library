clear all; close all; clc;

%% Input
fc      = 4;
fy      = 60;
H       = 20;
B       = 40;
cover   = 0.15*B;
rhosr   = 0.06;

nbB     = 2;
nbH     = 5;
axis    = 'y';
units   = 'US';

%% Define cross section object
conc_cross_section = Rectangle_Shape(H,B);

Ab = H*B*rhosr/(2*nbB+2*nbH-4);
reinforcement = reinf_rect(B-2*cover,H-2*cover,0,0,nbB,nbH,Ab);

section = RC(fc,fy,conc_cross_section,reinforcement,units);

%% Plot cross section
figure
section.plotSection
axis_limits('Margin',0.1)

%% Compute and plot interaction diagram
[P,M] = section.sectionInteraction2d(axis,'ACI');

figure
img = plot_data_on_image('data/Figure D.15.jpg');
img.set_point(1,[ 777 1637],[0.00 0.00])
img.set_point(2,[2321 1629],[0.55 0.00])
img.set_point(3,[ 768   92],[0.00 2.00])
img.show_image
img.plot_data(M/(fc*H*B*B),-P/(fc*H*B),'--xg','LineWidth',2);
