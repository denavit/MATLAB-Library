function [Ix,Iy] = find_intersection_between_two_lines(Ax,Ay,Bx,By,Cx,Cy,Dx,Dy)
% This function finds the intersection of lines AB and CD. The lines are 
% defined by two points but the intesection may occur and outside of the
% two points. If the lines are parallel, Ix and Iy are retured as empty
% matricies. 

if ( Ax == Bx && Cx == Dx )
    % Two vertical lines
    Ix = [];
    Iy = [];
elseif ( Ax == Bx )
    % AB is vertical
    CDm = (Dy-Cy)/(Dx-Cx);
    Ix = Ax;
    Iy = CDm*(Ix-Cx)+Cy;
elseif ( Cx == Dx )
    % CD is vertical
    ABm = (By-Ay)/(Bx-Ax);
    Ix = Cx;
    Iy = ABm*(Ix-Ax)+Ay;
else
    ABm = (By-Ay)/(Bx-Ax);
    CDm = (Dy-Cy)/(Dx-Cx);
    if ( ABm == CDm )
        % Two parallel lines
        Ix = [];
        Iy = [];
    else
        % Normal Situation
        Ix = (ABm*Ax-CDm*Cx-Ay+Cy)/(ABm-CDm);
        Iy = ABm*(Ix-Ax)+Ay;
    end
end
end