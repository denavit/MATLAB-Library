clear all; close all; clc;

load ShapeData_Rectangular_HSS.mat
numShapes = length(ShapeData_Rectangular_HSS);

label_aisc = {ShapeData_Rectangular_HSS(:).label};
H = [ShapeData_Rectangular_HSS(:).Ht];
B = [ShapeData_Rectangular_HSS(:).B];
t = [ShapeData_Rectangular_HSS(:).tdes];

A_aisc  = [ShapeData_Rectangular_HSS(:).A];
Ix_aisc = [ShapeData_Rectangular_HSS(:).Ix];
Iy_aisc = [ShapeData_Rectangular_HSS(:).Iy];
Sx_aisc = [ShapeData_Rectangular_HSS(:).Sx];
Sy_aisc = [ShapeData_Rectangular_HSS(:).Sy];
rx_aisc = [ShapeData_Rectangular_HSS(:).rx];
ry_aisc = [ShapeData_Rectangular_HSS(:).ry];
Zx_aisc = [ShapeData_Rectangular_HSS(:).Zx];
Zy_aisc = [ShapeData_Rectangular_HSS(:).Zy];
A_calc  = nan(1,numShapes);
Ix_calc = nan(1,numShapes);
Iy_calc = nan(1,numShapes);
Sx_calc = nan(1,numShapes);
Sy_calc = nan(1,numShapes);
rx_calc = nan(1,numShapes);
ry_calc = nan(1,numShapes);
Zx_calc = nan(1,numShapes);
Zy_calc = nan(1,numShapes);

for i = 1:numShapes
    shape = Rectangular_Tube_Shape(H(i),B(i),t(i),2*t(i));
    
    A_calc(i)  = shape.A;
    Ix_calc(i) = shape.I('strong');
    Iy_calc(i) = shape.I('weak');
    Sx_calc(i) = shape.S('strong');
    Sy_calc(i) = shape.S('weak');
    rx_calc(i) = shape.r('strong');
    ry_calc(i) = shape.r('weak');
    Zx_calc(i) = shape.Z('strong');
    Zy_calc(i) = shape.Z('weak');
end

%%
fprintf(' Property | Max Error |   Shape   | Calc Value | AISC Value \n')
fprintf('----------+-----------+-----------+------------+------------\n')
[err,ind] = max((A_calc-A_aisc)./A_aisc);
fprintf('%9s |%9.2f%% |%10s |%11g |%11g \n','A',100*err,label_aisc{ind},A_calc(ind),A_aisc(ind))
[err,ind] = max((Ix_calc-Ix_aisc)./Ix_aisc);
fprintf('%9s |%9.2f%% |%10s |%11g |%11g \n','Ix',100*err,label_aisc{ind},Ix_calc(ind),Ix_aisc(ind))
[err,ind] = max((Iy_calc-Iy_aisc)./Iy_aisc);
fprintf('%9s |%9.2f%% |%10s |%11g |%11g \n','Iy',100*err,label_aisc{ind},Iy_calc(ind),Iy_aisc(ind))
[err,ind] = max((Sx_calc-Sx_aisc)./Sx_aisc);
fprintf('%9s |%9.2f%% |%10s |%11g |%11g \n','Sx',100*err,label_aisc{ind},Sx_calc(ind),Sx_aisc(ind))
[err,ind] = max((Sy_calc-Sy_aisc)./Sy_aisc);
fprintf('%9s |%9.2f%% |%10s |%11g |%11g \n','Sy',100*err,label_aisc{ind},Sy_calc(ind),Sy_aisc(ind))
[err,ind] = max((rx_calc-rx_aisc)./rx_aisc);
fprintf('%9s |%9.2f%% |%10s |%11g |%11g \n','rx',100*err,label_aisc{ind},rx_calc(ind),rx_aisc(ind))
[err,ind] = max((ry_calc-ry_aisc)./ry_aisc);
fprintf('%9s |%9.2f%% |%10s |%11g |%11g \n','ry',100*err,label_aisc{ind},ry_calc(ind),ry_aisc(ind))
[err,ind] = max((Zx_calc-Zx_aisc)./Zx_aisc);
fprintf('%9s |%9.2f%% |%10s |%11g |%11g \n','Zx',100*err,label_aisc{ind},Zx_calc(ind),Zx_aisc(ind))
[err,ind] = max((Zy_calc-Zy_aisc)./Zy_aisc);
fprintf('%9s |%9.2f%% |%10s |%11g |%11g \n','Zy',100*err,label_aisc{ind},Zy_calc(ind),Zy_aisc(ind))
