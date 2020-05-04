clear all; close all; clc;

units = 'US';

% Concrete
H   = 60;
B   = 24;
fc  = 10;

% Steel 
shape_name = 'W36x256';
xs  = 0;
ys  = 6;
Fy  = 50;
Es  = 29000;

% Reinforcing
reinf_bar_size = '#10';
nbx = 4;
nby = 9;
dp  = 3;
Fylr = 60;

inlcude_compression_cap = false;

%% Setup

% Lookup steel data
shape_data = steel_shape_lookup(shape_name);
As = shape_data.A;
d  = shape_data.d;
tw = shape_data.tw;
bf = shape_data.bf;
tf = shape_data.tf;

% Lookup reinf data
reinf_bar_data = reinf_bar_lookup(reinf_bar_size);
Ab = reinf_bar_data.area;

% Reinf coordinates
num_bars = 2*(nbx+nby)-4;
yri = linspace(-H/2+dp,H/2-dp,nby);
yri = yri(2:(end-1));
xr = nan(num_bars,1);
yr = nan(num_bars,1);
xr(1:nbx) = linspace(-B/2+dp,B/2-dp,nbx);
yr(1:nbx) = -H/2+dp;
xr((nbx+1):(nbx+nby-2)) = B/2-dp;
yr((nbx+1):(nbx+nby-2)) = yri;
xr((nbx+nby-1):(2*nbx+nby-2)) = linspace(B/2-dp,-B/2+dp,nbx);
yr((nbx+nby-1):(2*nbx+nby-2)) = H/2-dp;
xr((2*nbx+nby-1):end) = -B/2+dp;
yr((2*nbx+nby-1):end) = fliplr(yri);

% Other items
Ec  = 57*sqrt(1000*fc);
Asr = Ab*num_bars;
Ac  = H*B - As - Asr;

id_steel = 1;
id_conc  = 2;
id_reinf = 3;

figsty = figureStyle('Display');

%% Create fiberSection object
fs = fiberSection;
fs.addMaterial(id_steel,Es,Fy);
fs.addMaterial( id_conc,Ec,0.85*fc);
fs.addMaterial(id_reinf,Es,Fylr);

% Steel Section
fs.addIShape(id_steel,d,tw,bf,tf,xs,ys)

% Concrete Section
b1 = tw/2;
b2 = bf/2;
b3 = B/2;
d1 = d/2-tf;
d2 = d/2;
d3 = H/2;
fs.addPatch('quad',id_conc,b3,d2+ys,-b3,d2+ys,-b3,d3,b3,d3);
fs.addPatch('quad',id_conc,b3,-d3,-b3,-d3,-b3,-d2+ys,b3,-d2+ys);
fs.addPatch('quad',id_conc,-b2+xs,d1+ys,-b3,d1+ys,-b3,d2+ys,-b2+xs,d2+ys);
fs.addPatch('quad',id_conc,b3,d1+ys,b2+xs,d1+ys,b2+xs,d2+ys,b3,d2+ys);
fs.addPatch('quad',id_conc,-b2+xs,-d2+ys,-b3,-d2+ys,-b3,-d1+ys,-b2+xs,-d1+ys);
fs.addPatch('quad',id_conc,b3,-d2+ys,b2+xs,-d2+ys,b2+xs,-d1+ys,b3,-d1+ys);
fs.addPatch('quad',id_conc,-b1+xs,-d1+ys,-b3,-d1+ys,-b3,d1+ys,-b1+xs,d1+ys);
fs.addPatch('quad',id_conc,b3,-d1+ys,b1+xs,-d1+ys,b1+xs,d1+ys,b3,d1+ys);

% Steel Reinforcement
for i = 1:num_bars
    fs.addFiber(id_reinf, Ab,xr(i),yr(i))
    fs.addFiber( id_conc,-Ab,xr(i),yr(i))
end

%% Create ACI_strain_compatibility object
scACI = ACI_strain_compatibility(fs);
%scACI.AxesOrigin = 'ScaledSectionCentroid';
%scACI.AxesOrigin = 'TransformedSectionCentroid';
%scACI.AxesOrigin = 'GrossSectionCentroid';

scACI.addConcreteBoundary(-B/2,-H/2,0);
scACI.addConcreteBoundary(-B/2, H/2,0);
scACI.addConcreteBoundary( B/2, H/2,0);
scACI.addConcreteBoundary( B/2,-H/2,0);
scACI.addSteelBoundary(-bf/2+xs,-d/2+ys,0);
scACI.addSteelBoundary(-bf/2+xs, d/2+ys,0);
scACI.addSteelBoundary( bf/2+xs, d/2+ys,0);
scACI.addSteelBoundary( bf/2+xs,-d/2+ys,0);
for i = 1:num_bars
    scACI.addSteelBoundary(xr(i),yr(i),0);
end

if inlcude_compression_cap
    scACI.maxCompressiveStrength = -0.85*(As*Fy + Asr*Fylr + 0.85*Ac*fc);
end

scACI.addMaterial(   'steel',id_steel,Fy,Es);
scACI.addMaterial('concrete',id_conc,fc,units);
scACI.addMaterial(   'steel',id_reinf,Fylr,Es);


%% Plot Section
figsty.picture(5,5);
fs.plotPatches;
fs.plotFibers;

%% Calculate and Plot Interaction Diagram
angle = 0;
numPoints = 50;
[P,Mx,My] = scACI.interactionSweep(angle,numPoints);

figsty.figure(5,5);
figsty.axes();
plot(Mx/12,-P)
plot(My/12,-P)
xlabel('Bending Moment (kip-ft)')
ylabel('Axial Load (kips, compression positive)')
legend('x-axis','y-axis')
