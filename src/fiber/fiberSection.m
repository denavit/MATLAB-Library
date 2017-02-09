classdef fiberSection < handle
   
    properties
        DefaultAxesOrigin = 'AsDefined';
    end
    
    properties (Access = private)
        data_m = [];
        data_A = [];
        data_z = [];
        data_y = [];
        dataBuilt = false;
        fibers_m = zeros(0,1);
        fibers_A = zeros(0,1);
        fibers_z = zeros(0,1);
        fibers_y = zeros(0,1);
        patches = cell(0,1);
        materials_id = zeros(0,1);
        materials_E  = zeros(0,1);
        materials_SF = zeros(0,1);
        nfy
        nfz
    end
    
    methods
        %% Constructor
        function obj = fiberSection(nfy,nfz)
            if nargin == 0
                obj.nfy = 100;
                obj.nfz = 100;
            else
                obj.nfy = nfy;
                obj.nfz = nfz;
            end
        end
        
        %% Change fiber discretization
        function setFiberDensity(obj,nfy,nfz)
            obj.nfy = nfy;
            obj.nfz = nfz;
            obj.dataBuilt = false;
        end         
        
        %% Define section - basic functions
        function addFiber(obj,matID,area,zLoc,yLoc)
            assert(isequal(size(matID),size(area),size(zLoc),size(yLoc)),...
                'matID, area, zLoc, and yLoc must all be the same size');
            assert(isnumeric(matID)&&isnumeric(area)&&isnumeric(zLoc)&&isnumeric(yLoc),...
                'matID, area, zLoc, and yLoc must be numeric');
            assert(iscolumn(matID),...
                'matID, area, zLoc, and yLoc must column vectors')
            obj.fibers_m = vertcat(obj.fibers_m,matID);
            obj.fibers_A = vertcat(obj.fibers_A,area);
            obj.fibers_z = vertcat(obj.fibers_z,zLoc);
            obj.fibers_y = vertcat(obj.fibers_y,yLoc);
            obj.dataBuilt = false;
        end     
        function addPatch(obj,type,varargin)
            switch lower(type)
                case 'quad'
                    newPatch = fiberPatch_quad(varargin{:});
                case 'circ'
                    newPatch = fiberPatch_circ(varargin{:});
                case 'fillet'
                    newPatch = fiberPatch_fillet(varargin{:});                    
                otherwise 
                    error('unknown patch type');
            end
            obj.patches = vertcat(obj.patches,{newPatch});
            obj.dataBuilt = false;
        end
        function addMaterial(obj,id,E,SF)
            if nargin < 3
                E = NaN;
            end
            if nargin < 4
                SF = NaN;
            end
            if isempty(E)
                E = NaN;
            end
            assert(isnumeric(id) && isscalar(id),'id should be a numeric scalar')
            assert(sum(obj.materials_id==id)==0,'Material Already Defined');
            obj.materials_id = vertcat(obj.materials_id,id);
            assert(isnumeric(E) && isscalar(E),'E should be a numeric scalar')
            obj.materials_E  = vertcat(obj.materials_E,E);
            assert(isnumeric(SF) && isscalar(SF),'SF should be a numeric scalar')
            obj.materials_SF = vertcat(obj.materials_SF,SF);
        end
        
        
        %% Information and Output
        function num = numFibers(obj)
            num = length(obj.fibers_m);
        end
        function num = numPatches(obj)
            num = length(obj.patches);
        end
        function mats = matIDs(obj)
            mat = obj.fibers_m;            
            for iPatch = 1:obj.numPatches
                mat = vertcat(mat,obj.patches{iPatch}.matIDs);
            end
            mats = unique(mat);
        end
        function [m,A,z,y] = fiberData(obj,AxesOrigin)
            if (obj.dataBuilt == false)
                obj.buildFiberData
            end
            if nargin < 2
                AxesOrigin = obj.DefaultAxesOrigin;
            end
            if strcmpi(AxesOrigin,'AsDefined')
                zo = 0;
                yo = 0;
            else
                [zo,yo] = obj.axesOrigin(AxesOrigin);
            end
            m = obj.data_m;
            A = obj.data_A;
            z = obj.data_z - zo;
            y = obj.data_y - yo;
        end
        function [z,y] = axesOrigin(obj,type)
            [m,A,z,y] = obj.fiberData('AsDefined');
            switch lower(type)
                case 'grosssectioncentroid'
                    sumA  = sum(A);
                    sumAz = sum(A.*z);
                    sumAy = sum(A.*y);
                    z = sumAz/sumA;
                    y = sumAy/sumA;
                case 'transformedsectioncentroid'
                    E = nan(size(m));
                    for i = 1:length(obj.materials_id)
                        E(m==obj.materials_id(i)) = obj.materials_E(i);
                    end
                    assert(~any(isnan(E)),'The stiffness of some materials is not defined')
                    sumEA  = sum(E.*A);
                    sumEAz = sum(E.*A.*z);
                    sumEAy = sum(E.*A.*y);
                    z = sumEAz/sumEA;
                    y = sumEAy/sumEA;
                case 'scaledsectioncentroid'
                    SF = nan(size(m));
                    for i = 1:length(obj.materials_id)
                        SF(m==obj.materials_id(i)) = obj.materials_SF(i);
                    end
                    assert(~any(isnan(SF)),'The scale factor of some materials is not defined')
                    sumSFA  = sum(SF.*A);
                    sumSFAz = sum(SF.*A.*z);
                    sumSFAy = sum(SF.*A.*y);
                    z = sumSFAz/sumSFA;
                    y = sumSFAy/sumSFA;
                otherwise
                    error('Unknown Axes Origin Type: %s',type)
            end            
        end
        function [zmin,zmax,ymin,ymax] = bounds(obj,AxesOrigin)
            if nargin < 2
                AxesOrigin = obj.DefaultAxesOrigin;
            end
            if strcmpi(AxesOrigin,'AsDefined')
                zo = 0;
                yo = 0;
            else
                [zo,yo] = obj.axesOrigin(AxesOrigin);
            end
            
            zmin =  Inf;
            zmax = -Inf;
            ymin =  Inf;
            ymax = -Inf;         
            
            if obj.numFibers > 0
                zmin = min(zmin,min(obj.fibers_z-zo));
                zmax = max(zmax,max(obj.fibers_z-zo));
                ymin = min(ymin,min(obj.fibers_y-yo));
                ymax = max(ymax,max(obj.fibers_y-yo));
            end
            
            for iPatch = 1:obj.numPatches
                [zminp,zmaxp,yminp,ymaxp] = obj.patches{iPatch}.bounds;
                zmin = min(zmin,zminp-zo);
                zmax = max(zmax,zmaxp-zo);
                ymin = min(ymin,yminp-yo);
                ymax = max(ymax,ymaxp-yo);
            end
        end
        function [ymax,ymin] = boundsAtAngle(obj,angle,AxesOrigin)
            if nargin < 2
                AxesOrigin = obj.DefaultAxesOrigin;
            end
            [~,~,z,y] = obj.fiberData(AxesOrigin);
            dist = [-sin(angle) cos(angle)]*[z y]';
            ymax = max(dist);
            ymin = min(dist);  
        end
        
        %% Print Function
        function printSectionProperties(obj,AxesOrigin,fid)
            if nargin < 2
                AxesOrigin = obj.DefaultAxesOrigin;
            end
            if nargin < 3
                fid = 1;
            end
            mats = obj.matIDs;
            [mat,A,z,y] = obj.fiberData(AxesOrigin);
            ATotal = 0;
            IzTotal = 0;
            IyTotal = 0;
            for i = 1:length(mats)
                ifibers = find(mat == mats(i));
                AMat  = sum(A(ifibers));
                IzMat = sum(A(ifibers).*y(ifibers).^2);
                IyMat = sum(A(ifibers).*z(ifibers).^2); 
                fprintf(fid,'Material ID %g:, Area = %g, Iz = %g, Iy = %g\n',...
                    mats(i),AMat,IzMat,IyMat);
                ATotal  = ATotal  + AMat;
                IzTotal = IzTotal + IzMat;
                IyTotal = IyTotal + IyMat;
            end
            fprintf(fid,'Total: Area = %g, Iz = %g, Iy = %g\n',ATotal,IzTotal,IyTotal);            
        end     
        
        %% Plotting Functions
        function plotSection(obj,AxesOrigin)
            if nargin < 2
                AxesOrigin = obj.DefaultAxesOrigin;
            end
            figure
            obj.plotFibers(AxesOrigin);
            obj.setAxisLimits(0.1);
            axis equal
            xlabel('Z Location')
            ylabel('Y Location')
        end
        function plotFibers(obj,AxesOrigin,id,color,size)
            if nargin < 2
                AxesOrigin = obj.DefaultAxesOrigin;
            end
            [m,~,z,y] = obj.fiberData(AxesOrigin);
            switch nargin
                case {1,2}
                    scatter(z,y,[],m);
                case 3
                    ind = m == id;
                    scatter(z(ind),y(ind),[],m(ind));
                otherwise
                    ind = m == id;
                    plot(z(ind),y(ind),'o','MarkerSize',size,'MarkerFaceColor',color,'Color',color);
            end
        end
        function plotPatches(obj,lineWidth)
            if nargin < 2
                lineWidth = 2;
            end
            hold all;
            for iPatch = 1:obj.numPatches
                obj.patches{iPatch}.plot(lineWidth);
            end
        end
        function setAxisLimits(obj,margin)
            [zmin,zmax,ymin,ymax] = obj.bounds;
            xlim([zmin-margin*(zmax-zmin) zmax+margin*(zmax-zmin)]);
            ylim([ymin-margin*(ymax-ymin) ymax+margin*(ymax-ymin)]);           
        end
              
        
        %% Define section - derived functions        
        function addStraightLayer(obj,matID,numFiber,areaFiber,zStart,yStart,zEnd,yEnd,skip)
            if nargin < 9
                skip = 'None';
            end
            matIDs = matID * ones(numFiber,1);
            areaFibers = areaFiber * ones(numFiber,1);
            zLocs = linspace(zStart,zEnd,numFiber)';
            yLocs = linspace(yStart,yEnd,numFiber)';
            switch lower(skip)
                case 'none'
                    obj.addFiber(matIDs,areaFibers,zLocs,yLocs);
                case 'i'
                    obj.addFiber(matIDs(2:end),areaFibers(2:end),zLocs(2:end),yLocs(2:end));
                case 'j'
                    obj.addFiber(matIDs(1:end-1),areaFibers(1:end-1),zLocs(1:end-1),yLocs(1:end-1));
                case 'both'
                    obj.addFiber(matIDs(2:end-1),areaFibers(2:end-1),zLocs(2:end-1),yLocs(2:end-1));
                otherwise
                    error('Unknown skip: %s',skip)
            end
            
        end
        function addCircLayer(obj,matID,numFiber,areaFiber,zCenter,yCenter,radius,startAng,endAng)
            if (nargin == 7)
                startAng = 0;
                endAng = 2*pi;
            end
            matIDs = matID * ones(numFiber,1);
            areaFibers = areaFiber * ones(numFiber,1);
            radii = radius * ones(numFiber,1);
            if mod(startAng,2*pi) == mod(endAng,2*pi)
                angles = linspace(startAng,endAng,numFiber+1)';
                angles = angles(1:numFiber);
            else
                angles = linspace(startAng,endAng,numFiber)';
            end
            [zLocs,yLocs] = pol2cart(angles,radii);
            zLocs = zLocs+zCenter;
            yLocs = yLocs+yCenter;
            obj.addFiber(matIDs,areaFibers,zLocs,yLocs);
        end          
        function addIShape(obj,matID,d,tw,bf,tf,zCenter,yCenter,angle,negative)
            if nargin < 7
                zCenter = 0;
                yCenter = 0;
            end
            if nargin < 9
                angle = 0;
            end
            if nargin < 10
                negative = false;
            end
            zPts = [ -bf/2 bf/2 ...
                -bf/2 -tw/2 tw/2 bf/2 ...
                -bf/2 -tw/2 tw/2 bf/2 ...
                    -bf/2 bf/2];
            yPts = [ d/2 d/2 ...
                d/2-tf  d/2-tf  d/2-tf  d/2-tf ...
               -d/2+tf -d/2+tf -d/2+tf -d/2+tf...
                    -d/2 -d/2];
            angle = deg2rad(angle);
            T = [ cos(angle) -sin(angle)
                  sin(angle)  cos(angle)];
            for i = 1:length(zPts);
                iPt = [zPts(i) yPts(i)]';
                iPtT = T*iPt;
                zPts(i) = iPtT(1);
                yPts(i) = iPtT(2);
            end
            zPts = zPts+zCenter;
            yPts = yPts+yCenter;
            obj.addPatch('quad',matID,zPts(1),yPts(1),zPts(2),yPts(2),zPts(6),yPts(6),zPts(3),yPts(3),negative);
            obj.addPatch('quad',matID,zPts(4),yPts(4),zPts(5),yPts(5),zPts(9),yPts(9),zPts(8),yPts(8),negative);
            obj.addPatch('quad',matID,zPts(7),yPts(7),zPts(10),yPts(10),zPts(12),yPts(12),zPts(11),yPts(11),negative);
        end
        function addWFShape(obj,matID,d,tw,bf,tf,k,negative)
            if nargin < 8
                negative = false;
            end
            zPts = [ -bf/2 bf/2 ...
                -bf/2 -tw/2 tw/2 bf/2 ...
                -bf/2 -tw/2 tw/2 bf/2 ...
                    -bf/2 bf/2];
            yPts = [ d/2 d/2 ...
                d/2-tf  d/2-tf  d/2-tf  d/2-tf ...
               -d/2+tf -d/2+tf -d/2+tf -d/2+tf...
                    -d/2 -d/2];
            obj.addPatch('quad',matID,zPts(1),yPts(1),zPts(2),yPts(2),zPts(6),yPts(6),zPts(3),yPts(3),negative);
            obj.addPatch('quad',matID,zPts(4),yPts(4),zPts(5),yPts(5),zPts(9),yPts(9),zPts(8),yPts(8),negative);
            obj.addPatch('quad',matID,zPts(7),yPts(7),zPts(10),yPts(10),zPts(12),yPts(12),zPts(11),yPts(11),negative);
            rad = k-tf;
            if rad > 0
                obj.addPatch('fillet',matID,zPts(4),yPts(4),rad,3,negative);
                obj.addPatch('fillet',matID,zPts(5),yPts(5),rad,4,negative);
                obj.addPatch('fillet',matID,zPts(8),yPts(8),rad,2,negative);
                obj.addPatch('fillet',matID,zPts(9),yPts(9),rad,1,negative);
            end
        end        
        function addRectHSS(obj,matID,D,B,t)
            ro = 2*t;
            ri = t;
            b1 = B/2-ro;
            b2 = B/2-(ro-ri);
            b3 = B/2;
            d1 = D/2-ro;
            d2 = D/2-(ro-ri);
            d3 = D/2;
            % Steel Section
            obj.addPatch('quad',matID,-b1,d2,-b1,d3,b1,d3,b1,d2);
            obj.addPatch('quad',matID,-b1,-d3,-b1,-d2,b1,-d2,b1,-d3);
            obj.addPatch('quad',matID,b2,-d1,b2,d1,b3,d1,b3,-d1);
            obj.addPatch('quad',matID,-b3,-d1,-b3,d1,-b2,d1,-b2,-d1);
            obj.addPatch('quad',matID,b1,d1,ri,ro,0.0,90.0);
            obj.addPatch('quad',matID,-b1,d1,ri,ro,90.0,180.0);
            obj.addPatch('quad',matID,-b1,-d1,ri,ro,180.0,270.0);
            obj.addPatch('quad',matID,b1,-d1,ri,ro,270.0,360.0);                     
        end        
        function addCCFT(obj,matID_steel,matID_conc,D,t)
            ro = D/2;
            ri = ro-t;
            obj.addPatch('circ',matID_steel,0.0,0.0,ri,ro);
            obj.addPatch('circ',matID_conc,0.0,0.0,0.0,ri);  
        end
        function addRCFT(obj,matID_steel,matID_conc,D,B,t,ri)
            if ri == 0
                b1 = B/2-t;
                b2 = B/2;
                d1 = D/2-t;
                d2 = D/2;
                % Steel Section
                obj.addPatch('quad',matID_steel,-b1, d1,-b1, d2, b1, d2, b1, d1);
                obj.addPatch('quad',matID_steel,-b1,-d2,-b1,-d1, b1,-d1, b1,-d2);
                obj.addPatch('quad',matID_steel,-b2,-d2,-b2, d2,-b1, d2,-b1,-d2);
                obj.addPatch('quad',matID_steel, b1,-d2, b1, d2, b2, d2, b2,-d2);
                % Concrete Section
                obj.addPatch('quad',matID_conc ,-b1,-d1,-b1, d1, b1, d1, b1,-d1);               
            else
                ro = ri+t;
                b1 = B/2-ro;
                b2 = B/2-(ro-ri);
                b3 = B/2;
                d1 = D/2-ro;
                d2 = D/2-(ro-ri);
                d3 = D/2;
                % Steel Section
                obj.addPatch('quad',matID_steel,-b1, d2,-b1, d3, b1, d3, b1, d2);
                obj.addPatch('quad',matID_steel,-b1,-d3,-b1,-d2, b1,-d2, b1,-d3);
                obj.addPatch('quad',matID_steel, b2,-d1, b2, d1, b3, d1, b3,-d1);
                obj.addPatch('quad',matID_steel,-b3,-d1,-b3, d1,-b2, d1,-b2,-d1);
                obj.addPatch('circ',matID_steel, b1, d1,ri,ro,   0.0,  pi/2);
                obj.addPatch('circ',matID_steel,-b1, d1,ri,ro,  pi/2,    pi);
                obj.addPatch('circ',matID_steel,-b1,-d1,ri,ro,    pi,3*pi/2);
                obj.addPatch('circ',matID_steel, b1,-d1,ri,ro,3*pi/2,  2*pi);
                % Concrete Section
                obj.addPatch('quad',matID_conc,-b2,-d1,-b2, d1, b2, d1, b2,-d1);
                obj.addPatch('quad',matID_conc,-b1, d1,-b1, d2, b1, d2, b1, d1);
                obj.addPatch('quad',matID_conc,-b1,-d2,-b1,-d1, b1,-d1, b1,-d2);
                obj.addPatch('circ',matID_conc, b1, d1,0.0,ri,   0.0,  pi/2);
                obj.addPatch('circ',matID_conc,-b1, d1,0.0,ri,  pi/2,    pi);
                obj.addPatch('circ',matID_conc,-b1,-d1,0.0,ri,    pi,3*pi/2);
                obj.addPatch('circ',matID_conc, b1,-d1,0.0,ri,3*pi/2,  2*pi);
            end
        end
        function addSRC(obj,matID_steel,matID_conc,matID_reinf,H,B,d,tw,bf,tf,zReinf,yReinf,dbReinf)
            % Steel Section
            obj.addIShape(matID_steel,d,tw,bf,tf)
            % Concrete Section
            b1 = tw/2;
            b2 = bf/2;
            b3 = B/2;
            d1 = d/2-tf;
            d2 = d/2;
            d3 = H/2;
            obj.addPatch('quad',matID_conc,b3,d2,-b3,d2,-b3,d3,b3,d3);
            obj.addPatch('quad',matID_conc,b3,-d3,-b3,-d3,-b3,-d2,b3,-d2);
            obj.addPatch('quad',matID_conc,-b2,d1,-b3,d1,-b3,d2,-b2,d2);
            obj.addPatch('quad',matID_conc,b3,d1,b2,d1,b2,d2,b3,d2);
            obj.addPatch('quad',matID_conc,-b2,-d2,-b3,-d2,-b3,-d1,-b2,-d1);
            obj.addPatch('quad',matID_conc,b3,-d2,b2,-d2,b2,-d1,b3,-d1);
            obj.addPatch('quad',matID_conc,-b1,-d1,-b3,-d1,-b3,d1,-b1,d1);
            obj.addPatch('quad',matID_conc,b3,-d1,b1,-d1,b1,d1,b3,d1);
            % Steel Reinforcement
            Ab = (pi/4)*dbReinf^2;
            for i = 1:length(zReinf)
                obj.addFiber(matID_reinf,Ab,zReinf(i),yReinf(i))
                obj.addFiber(matID_conc,-Ab,zReinf(i),yReinf(i))
            end
        end        
    end
    
    methods (Access = private)
        function buildFiberData(obj)
           
            % Find outer boundries of section
            [zMin,zMax,yMin,yMax] = obj.bounds;
            sfz = (zMax-zMin)/obj.nfz;
            sfy = (yMax-yMin)/obj.nfy;
            
            % Fibers
            obj.data_m = obj.fibers_m;
            obj.data_A = obj.fibers_A;
            obj.data_z = obj.fibers_z;
            obj.data_y = obj.fibers_y;

            % Patches 
            for iPatch = 1:obj.numPatches
                [m,A,z,y] = obj.patches{iPatch}.fiberData3d(sfz,sfy);
                obj.data_m = vertcat(obj.data_m,m);
                obj.data_A = vertcat(obj.data_A,A);
                obj.data_z = vertcat(obj.data_z,z);
                obj.data_y = vertcat(obj.data_y,y);
            end
            
            % Set dataBuilt flag
            obj.dataBuilt = true;
        end
    end
end

