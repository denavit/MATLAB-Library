classdef plastic_stress_distribution < handle
    
    properties
        AxesOrigin = 'AsDefined';
        fiberSection 
        materials_id = zeros(0,1);
        materials_Ft = zeros(0,1);
        materials_Fc = zeros(0,1);
    end

    methods
        %% Constructor
        function obj = plastic_stress_distribution(fiberSection,AxesOrigin)
            if nargin > 0
                obj.fiberSection = fiberSection;
            end
            if nargin > 1
                obj.AxesOrigin = AxesOrigin;
            end
        end
        
        %% Define section
        function addMaterial(obj,id,Ft,Fc)
            assert(sum(obj.materials_id==id)==0,'Material Already Defined');
            obj.materials_id = vertcat(obj.materials_id,id);
            obj.materials_Ft = vertcat(obj.materials_Ft,Ft);
            obj.materials_Fc = vertcat(obj.materials_Fc,Fc);
        end
        function set.fiberSection(obj,fiberSection)
            assert(isa(fiberSection,'fiberSection'),...
                'fiberSection should be a fiberSection object');
            obj.fiberSection = fiberSection;
        end
        
        %% Computation
        function [P,Mz,My] = computePoint(obj,zPoint,yPoint,angle)
            mats = obj.fiberSection.matIDs;
            [mat,A,z,y] = obj.fiberSection.fiberData(obj.AxesOrigin);
            
            a = sin(angle);
            b = -cos(angle);
            c = -sin(angle)*zPoint + cos(angle)*yPoint;
            P = 0;
            Mz = 0;
            My = 0;
            for i = 1:length(mats)
                ind = find(obj.materials_id==mats(i));
                assert(~isempty(ind),'Material %i is undefined',mats(i));
                Ft = obj.materials_Ft(ind);
                Fc = obj.materials_Fc(ind);
                tensFibers = find(mat(:,1) == mats(i) & a*z+b*y+c < 0);
                compFibers = find(mat(:,1) == mats(i) & a*z+b*y+c > 0);
                P  = P  + sum(Ft*A(tensFibers)) + sum(Fc*A(compFibers));
                Mz = Mz - sum(Ft*A(tensFibers).*y(tensFibers)) - ...
                    sum(Fc*A(compFibers).*y(compFibers));
                My = My + sum(Ft*A(tensFibers).*z(tensFibers)) + ...
                    sum(Fc*A(compFibers).*z(compFibers));
            end
        end      
        function [P,Mz,My] = interactionSweep(obj,angle,numPoints)           
            % Find points
            [~,~,z,y] = obj.fiberSection.fiberData(obj.AxesOrigin);
            dist = [-sin(angle) cos(angle)]*[z y]';
            maxDist = max(dist);
            minDist = min(dist);
            points = linspace(minDist,maxDist,numPoints);
            points(1) = points(1)-1.0e-6*(maxDist-minDist);
            points(end) = points(end)+1.0e-6*(maxDist-minDist);

            % Compute interaction
            P  = zeros(2*numPoints,1);
            Mz = zeros(2*numPoints,1);
            My = zeros(2*numPoints,1);
            for i = 1:numPoints
                zPoint = -sin(angle)*points(i);
                yPoint =  cos(angle)*points(i);
                                
                [iP,iMz,iMy] = computePoint(obj,zPoint,yPoint,angle);
                P(i)  = iP;
                Mz(i) = iMz;
                My(i) = iMy;
                
                [iP,iMz,iMy] = computePoint(obj,zPoint,yPoint,angle+pi);
                P(numPoints+i)  = iP;
                Mz(numPoints+i) = iMz;
                My(numPoints+i) = iMy;
            end
        end
        function [P,Mz,My] = interaction3d(obj,numPoints,numAngles)            
            P  = zeros(2*numPoints,numAngles);
            Mz = zeros(2*numPoints,numAngles);
            My = zeros(2*numPoints,numAngles);
            angles = linspace(0,pi,numAngles);
            
            for i = 1:numAngles
                [iP,iMz,iMy] = obj.interactionSweep(angles(i),numPoints);
                P(:,i) = iP;
                Mz(:,i) = iMz;
                My(:,i) = iMy;
            end
        end
        function results = interaction3d_2(obj,numP,numAngles,numPoints_calc,numAngles_calc,type)
            
            if nargin < 4
                numPoints_calc = 50;
            end
            if nargin < 5
                numAngles_calc = 2*numAngles;
            end
            if nargin < 6
                type = 'full';
            end
            
            switch lower(type)
                case 'full'
                    angles_calc = linspace(0,2*pi,numAngles_calc+1);
                case 'quadrant'
                    angles_calc = linspace(0,pi/2,numAngles_calc);
                otherwise
                    error('Unknown type: %s');
            end
            
            Mz_calc = nan(numP,numAngles_calc);
            My_calc = nan(numP,numAngles_calc);             
                        
            [~,~,z,y] = obj.fiberSection.fiberData(obj.AxesOrigin);
            for i = 1:numAngles_calc
                dist = [-sin(angles_calc(i)) cos(angles_calc(i))]*[z y]';
                maxDist     = max(dist);
                minDist     = min(dist);
                points      = linspace(minDist,maxDist,numPoints_calc);
                points(1)   = points(1)-1.0e-6*(maxDist-minDist);
                points(end) = points(end)+1.0e-6*(maxDist-minDist);
                
                iP  = nan(numPoints_calc,1);
                iMz = nan(numPoints_calc,1);
                iMy = nan(numPoints_calc,1);
                for j = 1:numPoints_calc
                    zPoint = -sin(angles_calc(i))*points(j);
                    yPoint =  cos(angles_calc(i))*points(j);
                    
                    [jP,jMz,jMy] = computePoint(obj,zPoint,yPoint,angles_calc(i)+pi);
                    iP(j)  = jP;
                    iMz(j) = jMz;
                    iMy(j) = jMy;
                end
                
                if i == 1
                    P_max  = iP(1);
                    Mz_max = iMz(1);
                    My_max = iMy(1);
                    P_min  = iP(end);
                    Mz_min = iMz(end);
                    My_min = iMy(end);
                    P = linspace(P_min,P_max,numP);
                end
                
                for j = 1:numP
                    [ind,x] = find_limit_point_in_vector(iP,P(j));
                    Mz_calc(j,i) = interpolate_vector(iMz,ind,x);
                    My_calc(j,i) = interpolate_vector(iMy,ind,x);
                end
            end

            results = struct;
            results.P  = nan(numP,numAngles);
            results.Mz = nan(numP,numAngles);
            results.My = nan(numP,numAngles);
            
            switch lower(type)
                case 'full'
                    angles = linspace(0,2*pi,numAngles+1);
                    angles = angles(1:(end-1));
                case 'quadrant'
                    angles = linspace(0,pi/2,numAngles);
                otherwise
                    error('Unknown type: %s');
            end            

            for i = 1:numP
                if i == 1
                    results.P(i,:)  = P_min;
                    results.Mz(i,:) = Mz_min;
                    results.My(i,:) = My_min;
                elseif i == numP
                    results.P(i,:)  = P_max;
                    results.Mz(i,:) = Mz_max;
                    results.My(i,:) = My_max;
                else
                    id = interactionDiagram2d(Mz_calc(i,:),My_calc(i,:));
                    d = id.radial_distance(angles);
                    results.P(i,:)  = P(i);
                    results.Mz(i,:) = d.*cos(angles);
                    results.My(i,:) = d.*sin(angles);
                end
            end
        end
    end
end

