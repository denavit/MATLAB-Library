classdef fiberPatch_circ < fiberPatch
    % obj=fiberPatch_circ(matID,zc,yc,ri,ro,a1,a2)
    %
    % Creates a set of fibers to describe a circular patch
    %
    % zCenter,yCenter - Coordinates of the center of the circle
    % intRad,extRad - inner and outer radius of the circular patch
    % nfRad - number of fibers in the radial direction
    % startAng,endAng - start and end angle of the circular patch
    %    (measured from the positive z axis)
    % nfArc - number of fibers in the circumferential direction

    properties
        matID
        zc
        yc
        ri
        ro
        a1 = 0;
        a2 = 2*pi;
        negative = false;
    end
    
    methods
        function obj=fiberPatch_circ(matID,zc,yc,ri,ro,varargin)
            obj.matID = matID;
            obj.zc = zc;
            obj.yc = yc;
            obj.ri = ri;
            obj.ro = ro;
            switch length(varargin)
                case 0
                    
                case 1
                    obj.negative = varargin{1};
                case 2
                    obj.a1 = varargin{1};
                    obj.a2 = varargin{2};
                case 3
                    obj.a1 = varargin{1};
                    obj.a2 = varargin{2};
                    obj.negative = varargin{3};
                otherwise
                    error('Bad number of arguments')
            end
        end
            
        function [m,A,y] = fiberData2d(obj,axis,sfy)
            assert(and(obj.a1==0,obj.a2==2*pi),'Function only implemented for a full circle')
            nf = ceil(obj.ro/sfy);
            dincr = obj.ro/nf;
            
            A = nan(nf,1);
            y = nan(nf,1);
            
            for i = 1:nf
                dfar  = obj.ro - (i-1)*dincr;
                dnear = obj.ro - i*dincr;
                if abs(dnear) < 1e-12*obj.ro 
                    dnear = 0;
                end
                if dnear >= obj.ri
                    [A1,y1] = circular_sector(dnear,obj.ro);
                    [A2,y2] = circular_sector( dfar,obj.ro);
                    nA = [A1 -A2];
                    ny = [y1 y2];
                elseif dfar >= obj.ri
                    [A1,y1] = circular_sector(dnear,obj.ro);
                    [A2,y2] = circular_sector( dfar,obj.ro);
                    [A3,y3] = circular_sector(dnear,obj.ri);
                    nA = [A1 -A2 -A3];
                    ny = [y1 y2 y3];
                else
                    [A1,y1] = circular_sector(dnear,obj.ro);
                    [A2,y2] = circular_sector( dfar,obj.ro);
                    [A3,y3] = circular_sector(dnear,obj.ri);
                    [A4,y4] = circular_sector( dfar,obj.ri);
                    nA = [A1 -A2 -A3 A4];
                    ny = [y1 y2 y3 y4];
                end
                A(i) = sum(nA);
                y(i) = sum(nA.*ny)/sum(nA);
            end

            % Get full circle
            A = vertcat(A, flipud(A));
            y = vertcat(y, flipud(y));
            
            % Adjust based on center of circle
            y = y + obj.yc;
            
            % Define material vector
            m = repmat(obj.matID,size(A));
            
            if obj.negative
                A = -A;
            end
        end
        function [m,A,z,y] = fiberData3d(obj,sfz,sfy)
            sfm = min([sfy sfz]);
                                    
            A = zeros(0,1);
            z = zeros(0,1);
            y = zeros(0,1);
            
            % Compute fiber data
            nfRad = max(1,ceil((obj.ro-obj.ri)/sfm));
            zRad = linspace(obj.ri,obj.ro,nfRad+1);
            
            for i = 1:nfRad
                iRi = zRad(i);
                iRo = zRad(i+1);
                nfArc = max(1,ceil(iRo*(obj.a2-obj.a1)/sfm));
                dAng = (obj.a2-obj.a1)/nfArc;
                zAng = linspace(obj.a1,obj.a2,nfArc+1); 
                iA = zeros(nfArc,1);
                iz = zeros(nfArc,1);
                iy = zeros(nfArc,1);
                for j = 1:nfArc
                    % Compute Area and Centroid
                    iArea1 =  (dAng/2)*iRo^2;
                    iArea2 = -(dAng/2)*iRi^2;
                    iArea = iArea1 + iArea2;
                    iY1 = 4*iRo*sin(dAng/2)/3/dAng;
                    iY2 = 4*iRi*sin(dAng/2)/3/dAng;
                    iY = (iArea1*iY1+iArea2*iY2)/iArea;
                    iAng = (zAng(j+1)+zAng(j))/2;
                    
                    % Assign fiber data                    
                    iA(j) = iArea;
                    iz(j) = obj.zc + iY*cos(iAng);
                    iy(j) = obj.yc + iY*sin(iAng);
                end
                A = vertcat(A,iA);
                z = vertcat(z,iz);
                y = vertcat(y,iy);
            end
            m = repmat(obj.matID,size(A));
            
            if obj.negative
                A = -A;
            end
        end 
        function [zmin,zmax,ymin,ymax] = bounds(obj)
            zmin = obj.zc-obj.ro;
            zmax = obj.zc+obj.ro;
            ymin = obj.yc-obj.ro;
            ymax = obj.yc+obj.ro;
        end
        function mat = matIDs(obj)
            mat = obj.matID;
        end
        function plot(obj,lineWidth)
            if nargin < 2
                lineWidth = 2;
            end
            if mod(obj.a2-obj.a1,2*pi) == 0
                angles = linspace(0,2*pi,100);
                hold all
                x = obj.zc + obj.ro*sin(angles);
                y = obj.yc + obj.ro*cos(angles);
                plot(x,y,'k-','LineWidth',lineWidth);
                if obj.ri ~= 0
                    x = obj.zc + obj.ri*sin(angles);
                    y = obj.yc + obj.ri*cos(angles);
                    plot(x,y,'k-','LineWidth',lineWidth);                    
                end
            else
                angles = obj.a1:pi/50:obj.a2;
                if obj.ri == 0
                    x = obj.zc + [0 obj.ro*sin(angles) 0];
                    y = obj.yc + [0 obj.ro*cos(angles) 0];
                    plot(x,y,'k-','LineWidth',lineWidth);                    
                else
                    x = obj.zc + [obj.ri*sin(angles) obj.ro*sin(fliplr(angles)) obj.ri*sin(obj.a1)];
                    y = obj.yc + [obj.ri*cos(angles) obj.ro*cos(fliplr(angles)) obj.ri*cos(obj.a1)];
                    plot(x,y,'k-','LineWidth',lineWidth);
                end
            end
            

        end
    end
end

function [A,y] = circular_sector(d,r)
if d < 0 || d > r
    error('d out of range');
end

if d == r
    A = 0;
    y = r;
else
    theta = 2*acos(d/r);
    A = 0.5*r^2*(theta-sin(theta));
    y = (4*r*sin(theta/2)^3)/(3*(theta-sin(theta)));
end
end