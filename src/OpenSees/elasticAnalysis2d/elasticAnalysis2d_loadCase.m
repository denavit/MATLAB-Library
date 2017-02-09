classdef elasticAnalysis2d_loadCase < OpenSeesAnalysis
    %elasticAnalysis2d_loadCase Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        model
        loadFactors
        
        
        constraint_type     = 'Transformation';
        constraint_alphaS   = [];
        constraint_alphaM   = [];
        numberer_type       = 'Plain';
        system_type         = 'UmfPack';
        test_type           = 'NormUnbalance';
        test_tol            = 1e-3;
        test_iter           = 10;
        test_flag           = 0;
        algorithm_type      = 'Newton';
        num_steps           = 10;
    end
    
    properties (SetAccess = private)
        results = struct;
    end
    
    methods
        %% Constructor
        function obj=elasticAnalysis2d_loadCase(model,loadFactors)
            % elasticAnalysis2d_loadCase(model,loadFactors,analysisOptions)
            obj.model = model;
            obj.loadFactors = loadFactors;
        end
        
        %% Set Functions
        function set.model(obj,model)
            if isa(model,'elasticAnalysis2d_model') && isscalar(model)
                obj.model = model;
            else
                error('model should be a elasticAnalysis2d_model');
            end
        end
        function set.loadFactors(obj,loadFactors)
            if isnumeric(loadFactors) && size(loadFactors,2)==2
                obj.loadFactors = loadFactors;
            else
                error('loadFactors should be numeric array with two columns');
            end
        end        
        
        %% Analysis Functions
        function analyze(obj,options)

            options.scratch_path = obj.scratchPath;
            options.loadFactors  = obj.loadFactors;
            
            % Filenames
            inputFilename = obj.scratchFile('elasticAnalysis2d_analysisDriver.tcl');
                       
            % Write Input Files
            fid = fopen(inputFilename,'w');
            
            % Source Model and Loading Input Files
            fprintf(fid, '# Define Model \n');
            modelFileData = obj.model.writeModelInputFile(fid,options);
            fprintf(fid, '\n');

            % Analysis Options
            fprintf(fid, '# Analysis Options \n');

            % Constraints
            switch obj.constraint_type
                case 'Plain'
                    fprintf(fid, 'constraints Plain \n');
                case 'Transformation'
                    fprintf(fid, 'constraints Transformation \n');
                case 'Penalty'
                    fprintf(fid, 'constraints Penalty %g %g \n',obj.constraint_alphaS,obj.constraint_alphaM);        
                case 'Lagrange'
                    fprintf(fid, 'constraints Lagrange %g %g \n',obj.constraint_alphaS,obj.constraint_alphaM);        
                otherwise
                    error('constraint_type not recgonized');
            end    

            % Numberer
            switch obj.numberer_type
                case 'Plain'
                    fprintf(fid, 'numberer Plain \n');
                case 'RCM'
                    fprintf(fid, 'numberer RCM \n');
                otherwise
                    error('numberer_type not recgonized');
            end   

            % System
            switch obj.system_type
                case 'BandGeneral'
                    fprintf(fid, 'system BandGeneral \n');
                case 'FullGeneral'
                    fprintf(fid, 'system FullGeneral \n');
                case 'UmfPack'
                    fprintf(fid, 'system UmfPack \n');
                otherwise
                    error('system_type not recgonized');
            end   

            % Test
            switch obj.test_type
                case 'NormUnbalance'
                    fprintf(fid, 'test NormUnbalance %g %i %i \n',obj.test_tol,obj.test_iter,obj.test_flag);
                case 'NormDispIncr'
                    fprintf(fid, 'test NormDispIncr %g %i %i \n',obj.test_tol,obj.test_iter,obj.test_flag);
                case 'EnergyIncr'
                    fprintf(fid, 'test EnergyIncr %g %i %i \n',obj.test_tol,obj.test_iter,obj.test_flag);
                otherwise
                    error('testType not recgonized');
            end   

            % Algorithm
            switch obj.algorithm_type
                case 'Newton'
                    fprintf(fid, 'algorithm Newton \n');
                case 'NewtonLineSearch'
                    fprintf(fid, 'algorithm NewtonLineSearch \n');
                case 'ModifiedNewton'
                    fprintf(fid, 'algorithm ModifiedNewton \n');
                otherwise
                    error('algorithm type not recgonized');
            end 
            fprintf(fid, '\n');

            % One step with no load
            fprintf(fid, '# Record Initial Conditions \n');
            fprintf(fid, 'record \n');
            fprintf(fid, '\n');

            % Apply load
            fprintf(fid, '# Apply Load \n');
            fprintf(fid, 'integrator LoadControl %g \n',1/obj.num_steps);
            fprintf(fid, 'analysis Static \n');
            fprintf(fid, 'set ok [analyze %i] \n',obj.num_steps);
            fprintf(fid, 'if {$ok == 0} { \n');
            fprintf(fid, '  # Analysis Successful \n');
            fprintf(fid, '  exit 1 \n');
            fprintf(fid, '} else { \n');
            fprintf(fid, '  # Analysis Failed \n');
            fprintf(fid, '  exit 2 \n');
            fprintf(fid, '} \n');
            fprintf(fid, '\n');

            % Close File
            fclose(fid);
            
            %% Run Analysis
            [status, result] = obj.runOpenSees(inputFilename);
            switch status
                case 1
                    % Analysis Successful
                case 2
                    fprintf('%s\n',result);
                    error('Analysis failed');
                otherwise
                    fprintf('%s\n',result);
                    error('Analysis ended in an unknown manner, exit code = %i',status);
            end
            
            % Read Results
            obj.results.textOutput = result;
            
            recorderFilenames = modelFileData.recorderFilenames;
            
            % Joint Displacements
            temp = dlmread(recorderFilenames.jointDisplacements);
            temp = temp(end,:);
            obj.results.jointDisplacements.x = temp(1:3:end);
            obj.results.jointDisplacements.y = temp(2:3:end);
            obj.results.jointDisplacements.r = temp(3:3:end);
            
            % Joint Reactions
            temp = dlmread(recorderFilenames.jointReactions);
            temp = temp(end,:);
            obj.results.jointReactions.x = temp(1:3:end);
            obj.results.jointReactions.y = temp(2:3:end);
            obj.results.jointReactions.r = temp(3:3:end);
            
            % Element Forces
            for i = 1:obj.model.numMembers
                numElements = obj.model.members(i).numElements;
                
                % Local forces for each member
                temp = dlmread(recorderFilenames.memberForces{i});
                temp = temp(end,:);
                temp2 = linspace(0,1,numElements+1);
                obj.results.memberForces(i).xi = zeros(1,2*numElements);
                obj.results.memberForces(i).xi(1:2:end) = temp2(1:end-1);
                obj.results.memberForces(i).xi(2:2:end) = temp2(2:end);
                obj.results.memberForces(i).P  = temp(1:3:end);
                obj.results.memberForces(i).P(1:2:end) = -obj.results.memberForces(i).P(1:2:end);
                obj.results.memberForces(i).V  = temp(2:3:end);
                obj.results.memberForces(i).V(2:2:end) = -obj.results.memberForces(i).V(2:2:end);
                obj.results.memberForces(i).M  = temp(3:3:end);
                obj.results.memberForces(i).M(1:2:end) = -obj.results.memberForces(i).M(1:2:end);
                
                % Displacements at the nodes for each member
                temp = dlmread(recorderFilenames.memberDisplacements{i});
                temp = temp(end,:);
                obj.results.memberDisplacements(i).xi = linspace(0,1,numElements+1);
                obj.results.memberDisplacements(i).x  = temp(1:3:end);
                obj.results.memberDisplacements(i).y  = temp(2:3:end);
                obj.results.memberDisplacements(i).r  = temp(3:3:end);
            end
            
            % Delete Scratch Files
            if obj.deleteFilesAfterAnalysis
                delete(inputFilename)
                delete(recorderFilenames.jointDisplacements);
                delete(recorderFilenames.jointReactions);
                for i = 1:obj.model.numMembers
                    delete(recorderFilenames.memberForces{i});
                    delete(recorderFilenames.memberDisplacements{i});
                end
            end
        end        
        
        %% Output Functions
        function memberForces = getMemberForces(obj,members)
            memberForces = struct([]);
            for i = 1:length(members);
                memberID = obj.model.findMember(members(i));
                if isempty(memberID)
                    memberForces(i).xi = [];
                    memberForces(i).P  = [];
                    memberForces(i).M  = [];
                    memberForces(i).V  = [];
                else
                    memberForces(i).xi = obj.results.memberForces(memberID).xi;
                    memberForces(i).P  = obj.results.memberForces(memberID).P;
                    memberForces(i).M  = obj.results.memberForces(memberID).M;
                    memberForces(i).V  = obj.results.memberForces(memberID).V;
                end
            end
        end
        function jointDisplacements = getJointDisplacements(obj,joints)
            jointDisplacements.x = obj.results.jointDisplacements.x(obj.model.findJoint(joints));
            jointDisplacements.y = obj.results.jointDisplacements.y(obj.model.findJoint(joints));
            jointDisplacements.r = obj.results.jointDisplacements.r(obj.model.findJoint(joints));     
        end
        function jointReactions = getJointReactions(obj,joints)
            jointReactions.x = obj.results.jointReactions.x(obj.model.findJoint(joints));
            jointReactions.y = obj.results.jointReactions.y(obj.model.findJoint(joints));
            jointReactions.r = obj.results.jointReactions.r(obj.model.findJoint(joints));            
        end     
        function disp = meanDisp(obj,storyJointIDs,direction)
            jointDisplacements = obj.getJointDisplacements(storyJointIDs);
            switch direction
                case 'x'
                    drifts = jointDisplacements.x;
                case 'y'
                    drifts = jointDisplacements.y;
                case 'r'
                    drifts = jointDisplacements.r;
                otherwise
                    error('Unknown direction');
            end
            disp = mean(drifts);
        end
        function [Fx,Fy,data] = forcesAtCut(obj,A,B,C)
			% @todo - should this function be based on initial coords or
			% displaced coords

            F = [0 0]';
            iMember = 1;
            for i = 1:obj.model.numMembers;
                memberID = obj.model.members(i).id;
                [iJointCoords,jJointCoords] = obj.model.memberCoords(memberID);
                iZ = A*iJointCoords(1) + B*iJointCoords(2) + C;
                jZ = A*jJointCoords(1) + B*jJointCoords(2) + C;
                if ( iZ == 0 || jZ == 0 )
                    warning('Cut is along a joint, empty matricies returned to avoid ambiguius results');
                    Fx = [];
                    Fy = [];
                    return;
                elseif ( iZ > 0 && jZ < 0 )
                    sideSign = 1;
                elseif ( iZ < 0 && jZ > 0 )
                    sideSign = -1;
                else
                    % Member is wholy on one side or the other
                    continue;
                end
                
                % Determine the location of the cut
                VEC = jJointCoords-iJointCoords;
                [theta,L] = cart2pol(VEC(1),VEC(2));
                if ( A == 0 )
                    [Ix,Iy] = find_intersection_between_two_lines(...
                        iJointCoords(1),iJointCoords(2),...
                        jJointCoords(1),jJointCoords(2),...
                        0.0,-C/B,1.0,-C/B);
                else
                    [Ix,Iy] = find_intersection_between_two_lines(...
                        iJointCoords(1),iJointCoords(2),...
                        jJointCoords(1),jJointCoords(2),...
                        -C/A,0.0,-(B+C)/A,1.0);
                end
                xi_cut = norm([Ix Iy]-iJointCoords)/L;
                
                % Determine axial load and shear at the cut
                memberForces = obj.getMemberForces(memberID);                
                axialLoadDiagram = piecewiseLinear(memberForces.xi,memberForces.P);
                P = axialLoadDiagram.val(xi_cut);
                shearDiagram = piecewiseLinear(memberForces.xi,memberForces.V);
                V = shearDiagram.val(xi_cut);
                
                % Transform to global coordinates
                T = [cos(theta) -sin(theta); sin(theta) cos(theta)];
                Fi = T*vertcat(sideSign*P,sideSign*-V);
                F = F + Fi;
                
                data.memberID(iMember) = memberID;
                data.Fx(iMember) = Fi(1);
                data.Fy(iMember) = Fi(2);
                iMember = iMember+1;
            end
            Fx = F(1);
            Fy = F(2);
            if nargout < 3
                clear data
            end
        end
        function print(obj,fid,flag)
            % print()
            % print(fid)
            % print(fid,flag)
            if (nargin < 2); fid = 1; end
            if (nargin < 3); flag = 1; end

            fprintf(fid,'elasticAnalysis2d Load Case:\n');
            fprintf(fid,'  Load Factors:\n');
            for i = 1:size(obj.loadFactors,1)
                fprintf(fid,'    Load Pattern = %i, Factor = %g\n',...
                    obj.loadFactors(i,1),obj.loadFactors(i,2));
            end
        end 
        
        %% Plotting Functions
        function plotDeformedJoints(obj,amp)
            for i = 1:length(obj.model.joints)
                plot(obj.model.joints(i).coords(1)+amp*obj.results.jointDisplacements.x(i),...
                    obj.model.joints(i).coords(2)+amp*obj.results.jointDisplacements.y(i),...
                    'ro','MarkerFaceColor','r','MarkerSize',7)
            end    
        end
        function plotDeformedMembers(obj,amp)
            for i = 1:length(obj.model.members)
                [iJointCoords,jJointCoords] = obj.model.memberCoords(i);
                memberID = obj.model.findMember(i);
                numElements = obj.model.members(memberID).numElements;
                xCoords = linspace(iJointCoords(1),jJointCoords(1),numElements+1);
                yCoords = linspace(iJointCoords(2),jJointCoords(2),numElements+1);
                xCoords = xCoords + amp*obj.results.memberDisplacements(i).x;
                yCoords = yCoords + amp*obj.results.memberDisplacements(i).y; 
                plot(xCoords,yCoords,...
                    'k-','LineWidth',2)
            end   
        end
        function plotBendingMoment(obj,amp)
            for i = 1:length(obj.model.members)          
                iJoint = obj.model.findJoint(obj.model.members(i).jointIDs(1));
                jJoint = obj.model.findJoint(obj.model.members(i).jointIDs(2));
                iJointCoords = obj.model.joints(iJoint).coords;
                jJointCoords = obj.model.joints(jJoint).coords; 
                VEC = jJointCoords-iJointCoords;
                [theta,L] = cart2pol(VEC(1),VEC(2));                
                X = [0 obj.results.memberForces(i).xi 1]*L;
                Y = [0 obj.results.memberForces(i).M 0]*amp;
                T = [cos(theta) -sin(theta); sin(theta) cos(theta)];
                XY = T*vertcat(X,Y);
                fill(iJointCoords(1)+XY(1,:),iJointCoords(2)+XY(2,:),[0 0 1]);
            end
        end
        function plotShear(obj,amp)
            for i = 1:length(obj.model.members)
                iJoint = obj.model.findJoint(obj.model.members(i).jointIDs(1));
                jJoint = obj.model.findJoint(obj.model.members(i).jointIDs(2));
                iJointCoords = obj.model.joints(iJoint).coords;
                jJointCoords = obj.model.joints(jJoint).coords;                
                VEC = jJointCoords-iJointCoords;
                [Q,L] = cart2pol(VEC(1),VEC(2));
                X = [0 obj.results.memberForces(i).xi 1]*L;
                Y = [0 obj.results.memberForces(i).V 0]*amp;
                T = [cos(Q) -sin(Q); sin(Q) cos(Q)];
                XY = T*vertcat(X,Y);
                fill(iJointCoords(1)+XY(1,:),iJointCoords(2)+XY(2,:),[0 0 1]);
            end
        end  
        function plotAxialLoad(obj,amp)
            for i = 1:length(obj.model.members)
                iJoint = obj.model.findJoint(obj.model.members(i).jointIDs(1));
                jJoint = obj.model.findJoint(obj.model.members(i).jointIDs(2));
                iJointCoords = obj.model.joints(iJoint).coords;
                jJointCoords = obj.model.joints(jJoint).coords; 
                VEC = jJointCoords-iJointCoords;
                [Q,L] = cart2pol(VEC(1),VEC(2));
                X = [0 obj.results.memberForces(i).xi 1]*L;
                Y = [0 obj.results.memberForces(i).P 0]*amp;
                T = [cos(Q) -sin(Q); sin(Q) cos(Q)];
                XY = T*vertcat(X,Y);
                fill(iJointCoords(1)+XY(1,:),iJointCoords(2)+XY(2,:),[0 0 1]);
            end
        end  
        
        function plotMemberLoads(obj,memberID)
            memberForces = obj.getMemberForces(memberID);
            figure
            subplot(3,1,1)
            plot(memberForces.xi,memberForces.P)
            ylabel('Axial Load')
            subplot(3,1,2)
            plot(memberForces.xi,memberForces.V)
            ylabel('Shear Force')
            subplot(3,1,3)
            plot(memberForces.xi,memberForces.M)
            ylabel('Bending Moment')
        end
        
        function plotStructure(obj,amp)
            figure
            hold all
            obj.model.plotJoints;
            obj.model.plotMembers;
            obj.plotDeformedJoints(amp);
            obj.plotDeformedMembers(amp);            
            %[xMin,xMax,yMin,yMax] = obj.model.structureBounds;
            %margin = 0.1;
            %xlim([xMin-margin*(xMax-xMin) xMax+margin*(xMax-xMin)])
            %ylim([yMin-margin*(yMax-yMin) yMax+margin*(yMax-yMin)])
        end
  
    end 
end

