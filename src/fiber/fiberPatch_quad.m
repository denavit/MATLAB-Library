classdef fiberPatch_quad < fiberPatch
    % obj=fiberPatch_quad(matID,vertexCoords)
    % obj=fiberPatch_quad(matID,zI,yI,zJ,yJ,zK,yK,zL,yL)
    %
    % Creates a set of fibers to describe a quadrilateral patch
    %
    % vertexCoords is a 4x2 matrix with the coodinates of the
    % verticies of the quadrilateral.
    %
    % - Clockwise ordering results in a positive area
    % - Counter-clockwise ordering results in a negative area
    %
    %   y |
    %     |    J  o-----o K
    %     |      /     /
    %     |     /     /
    %     |  I o-----o L
    %     |______________
    %                   z
    %
    
    properties
        matID
        zI
        yI
        zJ
        yJ
        zK
        yK
        zL
        yL
        negative = false;
    end
    
    methods
        function obj=fiberPatch_quad(matID,zI,yI,zJ,yJ,zK,yK,zL,yL,negative)
            obj.matID = matID;
            if nargin == 2 || nargin == 3
                assert(~isequal([4 2],size(zI)),'vertexCoords should be a 4x2 matrix');
                obj.zI = zI(1,1);
                obj.yI = yI(2,1);
                obj.zJ = zJ(1,2);
                obj.yJ = yJ(2,2);
                obj.zK = zK(1,3);
                obj.yK = yK(2,3);
                obj.zL = zL(1,4);
                obj.yL = yL(2,4);
                if nargin == 3
                    obj.negative = negative;
                end
            elseif nargin == 9 || nargin == 10
                obj.zI = zI;
                obj.yI = yI;
                obj.zJ = zJ;
                obj.yJ = yJ;
                obj.zK = zK;
                obj.yK = yK;
                obj.zL = zL;
                obj.yL = yL;
                if nargin == 10
                    obj.negative = negative;
                end
            else
                error('Bad number of arguments')
            end            
        end
            
        function [m,A,z] = fiberData2d(obj,axis,sfz)
            error('Not yet implemented');
        end
        function [m,A,z,y] = fiberData3d(obj,sfz,sfy)
            
            % Determine local fiber discretization
            [nfIJ,nfJK] = obj.localFiberDiscretization3d(sfz,sfy);
            
            patchCoords = obj.vertexCoords;
            
            m = repmat(obj.matID,nfIJ*nfJK,1);
            A = zeros(nfIJ*nfJK,1);
            z = zeros(nfIJ*nfJK,1);
            y = zeros(nfIJ*nfJK,1);
            
            % Compute fiber data
            deltaXi = 2.0/nfIJ;
            deltaEta = 2.0/nfJK;
            for j = 1:nfIJ
                for k = 1:nfJK
                    % Compute Cell Vertex Coordinates
                    cellVertexCoordsN = [
                        -1.0+deltaXi*(j-1) -1.0+deltaEta*(k-1)
                        -1.0+deltaXi*j     -1.0+deltaEta*(k-1)
                        -1.0+deltaXi*j     -1.0+deltaEta*k
                        -1.0+deltaXi*(j-1) -1.0+deltaEta*k ];
                    cellVertexCoords = zeros(4,2);
                    
                    for i = 1:4
                        xi  = cellVertexCoordsN(i,1);
                        eta = cellVertexCoordsN(i,2);
                        N = [(1.0-xi)*(1.0-eta) ...
                            (1.0+xi)*(1.0-eta) ...
                            (1.0+xi)*(1.0+eta) ...
                            (1.0-xi)*(1.0+eta)]/4.0;
                        cellVertexCoords(i,:) = N*patchCoords;
                    end
                    
                    % Compute Area and Centroid
                    l = (j-1)*nfJK + k;
                    [A(l),z(l),y(l)] = quadCell(cellVertexCoords);
                end
            end
            
            if obj.negative
                A = -A;
            end
        end
        function [zmin,zmax,ymin,ymax] = bounds(obj)
            zmin = min([obj.zI obj.zJ obj.zK obj.zL]);
            zmax = max([obj.zI obj.zJ obj.zK obj.zL]);
            ymin = min([obj.yI obj.yJ obj.yK obj.yL]);
            ymax = max([obj.yI obj.yJ obj.yK obj.yL]);
        end
        function mat = matIDs(obj)
            mat = obj.matID;
        end
        function plot(obj,lineWidth)
            if nargin < 2
                lineWidth = 2;
            end
            x = [obj.zI obj.zJ obj.zK obj.zL obj.zI];
            y = [obj.yI obj.yJ obj.yK obj.yL obj.yI];
            plot(x,y,'k-','LineWidth',lineWidth);
        end
        function [nfIJ,nfJK] = localFiberDiscretization3d(obj,sfy,sfz)
            nfIJ = ceil(max([1 ...
                abs(obj.zI-obj.zJ)/sfz abs(obj.yI-obj.yJ)/sfy ...
                abs(obj.zK-obj.zL)/sfz abs(obj.yK-obj.yL)/sfy]));
            nfJK = ceil(max([1 ...
                abs(obj.zJ-obj.zK)/sfz abs(obj.yJ-obj.yK)/sfy ...
                abs(obj.zL-obj.zI)/sfz abs(obj.yL-obj.yI)/sfy]));
        end
        function coords = vertexCoords(obj)
            coords = [
                obj.zI obj.yI 
                obj.zJ obj.yJ 
                obj.zK obj.yK
                obj.zL obj.yL ];
        end
    end
end