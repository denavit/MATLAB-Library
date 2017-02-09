clear all; close all; clc;

load ShapeData_Wide_Flange.mat
numShapes = length(ShapeData_Wide_Flange);

label_aisc = {ShapeData_Wide_Flange(:).label};
d  = [ShapeData_Wide_Flange(:).d];
tw = [ShapeData_Wide_Flange(:).tw];
bf = [ShapeData_Wide_Flange(:).bf];
tf = [ShapeData_Wide_Flange(:).tf];
k  = [ShapeData_Wide_Flange(:).kdes];

A_aisc  = [ShapeData_Wide_Flange(:).A];
Ix_aisc = [ShapeData_Wide_Flange(:).Ix];
Iy_aisc = [ShapeData_Wide_Flange(:).Iy];
Sx_aisc = [ShapeData_Wide_Flange(:).Sx];
Sy_aisc = [ShapeData_Wide_Flange(:).Sy];
rx_aisc = [ShapeData_Wide_Flange(:).rx];
ry_aisc = [ShapeData_Wide_Flange(:).ry];
Zx_aisc = [ShapeData_Wide_Flange(:).Zx];
Zy_aisc = [ShapeData_Wide_Flange(:).Zy];
J_aisc  = [ShapeData_Wide_Flange(:).J];
A_calc  = nan(1,numShapes);
Ix_calc = nan(1,numShapes);
Iy_calc = nan(1,numShapes);
Sx_calc = nan(1,numShapes);
Sy_calc = nan(1,numShapes);
rx_calc = nan(1,numShapes);
ry_calc = nan(1,numShapes);
Zx_calc = nan(1,numShapes);
Zy_calc = nan(1,numShapes);
J_calc  = nan(1,numShapes);

for i = 1:numShapes
    shape = I_Shape(d(i),tw(i),bf(i),tf(i),k(i));
    
    A_calc(i)  = shape.A;
    Ix_calc(i) = shape.I('strong');
    Iy_calc(i) = shape.I('weak');
    Sx_calc(i) = shape.S('strong');
    Sy_calc(i) = shape.S('weak');
    rx_calc(i) = shape.r('strong');
    ry_calc(i) = shape.r('weak');
    Zx_calc(i) = shape.Z('strong');
    Zy_calc(i) = shape.Z('weak');
    J_calc(i)  = shape.J;
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
[err,ind] = max((J_calc-J_aisc)./J_aisc);
fprintf('%9s |%9.2f%% |%10s |%11g |%11g \n','J',100*err,label_aisc{ind},J_calc(ind),J_aisc(ind))
