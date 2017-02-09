classdef fiberPatch_fillet < fiberPatch
    % obj=fiberPatch_fillet(matID,zv,yv,r,quadrant)
    %
    % Creates a set of fibers to describe a circular fillet
    %

    properties
        matID
        zv
        yv
        r
        quadrant
        negative = false;
    end
    
    methods
        function obj=fiberPatch_fillet(matID,zv,yv,r,quadrant,negative)
            obj.matID = matID;
            obj.zv = zv;
            obj.yv = yv;
            obj.r = r;
            obj.quadrant = quadrant;
            if nargin > 5
                obj.negative = negative;
            end
        end
            
        function [m,A,z] = fiberData2d(obj,axis,sfz)
            error('Not yet implemented');
        end
        function [m,A,z,y] = fiberData3d(obj,sfz,sfy)
            sfm = min([sfy sfz]);
            nf = max(1,ceil(obj.r/sfm));
            
            if (nf == 1)
                A = (1-pi/4)*obj.r^2;
                iCentroid = 2/(12-3*pi)*obj.r;
                switch obj.quadrant
                    case 1
                        z = obj.zv+iCentroid;
                        y = obj.yv+iCentroid;
                    case 2
                        z = obj.zv-iCentroid;
                        y = obj.yv+iCentroid;
                    case 3
                        z = obj.zv-iCentroid;
                        y = obj.yv-iCentroid;
                    case 4
                        z = obj.zv+iCentroid;
                        y = obj.yv-iCentroid;
                    otherwise
                        error('Unknown quadrant');
                end
            else
                switch obj.quadrant
                    case 1
                        zPtsA = [obj.zv*ones(1,nf) linspace(obj.zv,obj.zv+obj.r,nf+1)];
                        yPtsA = [linspace(obj.yv+obj.r,obj.yv,nf+1) obj.yv*ones(1,nf)];
                        angles = linspace(pi,3*pi/2,2*nf+1);
                        zPtsB = obj.zv+obj.r + obj.r*cos(angles);
                        yPtsB = obj.yv+obj.r + obj.r*sin(angles);
                    case 2
                        zPtsA = [linspace(obj.zv-obj.r,obj.zv,nf+1) obj.zv*ones(1,nf)];
                        yPtsA = [obj.yv*ones(1,nf) linspace(obj.yv,obj.yv+obj.r,nf+1)];
                        angles = linspace(3*pi/2,2*pi,2*nf+1);
                        zPtsB = obj.zv-obj.r + obj.r*cos(angles);
                        yPtsB = obj.yv+obj.r + obj.r*sin(angles);
                    case 3
                        zPtsA = [obj.zv*ones(1,nf) linspace(obj.zv,obj.zv-obj.r,nf+1)];
                        yPtsA = [linspace(obj.yv-obj.r,obj.yv,nf+1) obj.yv*ones(1,nf)];
                        angles = linspace(0,pi/2,2*nf+1);
                        zPtsB = obj.zv-obj.r + obj.r*cos(angles);
                        yPtsB = obj.yv-obj.r + obj.r*sin(angles);
                    case 4
                        zPtsA = [linspace(obj.zv+obj.r,obj.zv,nf+1) obj.zv*ones(1,nf)];
                        yPtsA = [obj.yv*ones(1,nf) linspace(obj.yv,obj.yv-obj.r,nf+1)];
                        angles = linspace(pi/2,pi,2*nf+1);
                        zPtsB = obj.zv+obj.r + obj.r*cos(angles);
                        yPtsB = obj.yv-obj.r + obj.r*sin(angles);
                    otherwise
                        error('Unknown quadrant');
                end
                
                A = zeros(0,1);
                z = zeros(0,1);
                y = zeros(0,1);
                
                for i = 1:(length(zPtsA)-1)
                    length1 = sqrt((zPtsA(i)-zPtsB(i))^2+(yPtsA(i)-yPtsB(i))^2);
                    length2 = sqrt((zPtsA(i+1)-zPtsB(i+1))^2+(yPtsA(i+1)-yPtsB(i+1))^2);
                    nf2 = ceil(max([1 max([length1 length2])/sfm]));
                    zLine1 = linspace(zPtsA(i),zPtsB(i),nf2+1);
                    yLine1 = linspace(yPtsA(i),yPtsB(i),nf2+1);
                    zLine2 = linspace(zPtsA(i+1),zPtsB(i+1),nf2+1);
                    yLine2 = linspace(yPtsA(i+1),yPtsB(i+1),nf2+1);
                    for j = 1:nf2
                        cellVertexCoords = [ zLine1(j)   yLine1(j)   ;
                            zLine1(j+1) yLine1(j+1) ;
                            zLine2(j+1) yLine2(j+1) ;
                            zLine2(j)   yLine2(j)   ];
                        [iA,iz,iy] = quadCell(cellVertexCoords);
                        A = vertcat(A,iA);
                        z = vertcat(z,iz);
                        y = vertcat(y,iy);
                    end
                end
            end
            m = repmat(obj.matID,size(A));
            
            if obj.negative
                A = -A;
            end
        end
        function [zmin,zmax,ymin,ymax] = bounds(obj)
            switch obj.quadrant
                case 1
                    zmin = obj.zv;
                    zmax = obj.zv+obj.r;
                    ymin = obj.yv;
                    ymax = obj.yv+obj.r;
                case 2
                    zmin = obj.zv-obj.r;
                    zmax = obj.zv;
                    ymin = obj.yv;
                    ymax = obj.yv+obj.r;                    
                case 3
                    zmin = obj.zv-obj.r;
                    zmax = obj.zv;
                    ymin = obj.yv-obj.r;
                    ymax = obj.yv;                    
                case 4
                    zmin = obj.zv;
                    zmax = obj.zv+obj.r;
                    ymin = obj.yv-obj.r;
                    ymax = obj.yv;
                otherwise
                    error('Unknown quadrant');
            end
        end
        function mat = matIDs(obj)
            mat = obj.matID;
        end
        function plot(obj,lineWidth)
            if nargin < 2
                lineWidth = 2;
            end
            switch obj.quadrant
                case 1
                    angles = linspace(pi,1.5*pi,50);
                    x = obj.zv + [0 obj.r*( 1+cos(angles)) 0];
                    y = obj.yv + [0 obj.r*( 1+sin(angles)) 0];
                case 2
                    angles = linspace(1.5*pi,2*pi,50);
                    x = obj.zv + [0 obj.r*(-1+cos(angles)) 0];
                    y = obj.yv + [0 obj.r*( 1+sin(angles)) 0];
                case 3
                    angles = linspace(0,0.5*pi,50);
                    x = obj.zv + [0 obj.r*(-1+cos(angles)) 0];
                    y = obj.yv + [0 obj.r*(-1+sin(angles)) 0];
                case 4
                    angles = linspace(0.5*pi,pi,50);
                    x = obj.zv + [0 obj.r*( 1+cos(angles)) 0];
                    y = obj.yv + [0 obj.r*(-1+sin(angles)) 0];
                otherwise
                    error('Unknown quadrant');
            end
            plot(x,y,'k-','LineWidth',lineWidth);
        end
    end
end

