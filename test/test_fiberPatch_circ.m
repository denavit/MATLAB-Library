clear all; close all; clc;

D = 12;
t = 0.25;

sfz = 0.5;
sfy = 0.5;

% Exact Formulas
sh = Circular_Tube_Shape(D,t);

A_exact  = sh.A;
Ix_exact = sh.I('srtong');

% Fiber Data
fp = fiberPatch_circ(1,0,0,D/2-t,D/2);
[m,A,z,y] = fp.fiberData3d(sfz,sfy);
%[m,A,y] = fp.fiberData2d('strong',sfy);

A_fiber  = sum(A);
Ix_fiber = sum(A.*y.^2);

%% Print Summary
fprintf('Area:\n');
fprintf('exact: %.3f \n',A_exact);
fprintf('fiber: %.3f \n',A_fiber);
fprintf('%%diff: %.2f%% \n',100*(A_fiber-A_exact)/A_exact);
fprintf('Moment of Intertia:\n');
fprintf('exact: %.3f \n',Ix_exact);
fprintf('fiber: %.3f \n',Ix_fiber);
fprintf('%%diff: %.2f%% \n',100*(Ix_fiber-Ix_exact)/Ix_exact);

%% Plot Fiber Data
figure
scatter(z,y)
axis equal