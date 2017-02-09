classdef mdofShearBuilding2d < OpenSeesAnalysis
    
    properties
        nStories
        
        storyMass
        storyStiffness
        storySpringDefinition
        
        damping_ModeA  = 1;
        damping_ModeB  = 3;
        damping_RatioA = 0.02;
        damping_RatioB = 0.02;
        
        controlStory = 'roof';

        pushover_stepSize   = 0.001;
        pushover_maxDrift   = 6.0;
    end
    
    methods
        %% Constructor
        function obj = mdofShearBuilding2d(nStories)
            % Constructor
            %
            % obj = mdofShearBuilding2d(nStories)
            %
            obj.nStories = nStories;
        end
        
        %% Set functions
        function set.storyMass(obj,storyMass)
            if ~isnumeric(storyMass)
                error('storyMass should be numeric');
            end
            if isvectorsize(storyMass,obj.nStories)
                obj.storyMass = storyMass;
            else 
                error('storyMass should be vector of length %i (number of stories)',obj.nStories);
            end
        end
        function set.storyStiffness(obj,storyStiffness)
            if ~isnumeric(storyStiffness)
                error('storyStiffness should be numeric');
            end
            if isvectorsize(storyStiffness,obj.nStories)
                obj.storyStiffness = storyStiffness;
            else 
                error('storyStiffness should be vector of length %i (number of stories)',obj.nStories);
            end
        end     
        function set.storySpringDefinition(obj,storySpringDefinition)
            if ~iscell(storySpringDefinition)
                error('storySpringDefinition should be a cell vector');
            end
            if isvectorsize(storySpringDefinition,obj.nStories)
                obj.storySpringDefinition = storySpringDefinition;
            else 
                error('storySpringDefinition should be cell vector of length %i (number of stories)',obj.nStories);
            end
        end          
        
        %% Functions
        function s = controlStory1(obj)
            if isequal(obj.controlStory,'roof')
                s = obj.nStories;
            else
                s = obj.controlStory;
            end
        end
        
        %% Analyses 
        function results = pushover(obj,F,type,varargin)
            assert(isnumeric(F) & isvectorsize(F,obj.nStories),...
                'F should be a numeric verctor of length %i (number of stories)',obj.nStories);
            if iscolumn(F)
                F = F';
            end

            % Initilize Results
            results = struct;
            results.F = F;            
            
            switch lower(type)
                case 'targetdrift'
                    targetDrift      = varargin{1};
                    
                    results.targetDrift = targetDrift;
                case 'targetpostpeakratio'
                    targetPostPeakRatio = varargin{1};
                    
                    results.targetPostPeakRatio = targetPostPeakRatio;
                otherwise
                    error('Unknown analysis type: %s',type);
            end
            
            % Filenames
            filename_input          = obj.scratchFile('mdofShearBuilding2d_input.tcl');
            filename_output_def     = obj.scratchFile('mdofShearBuilding2d_disp.out');
            filename_output_force   = obj.scratchFile('mdofShearBuilding2d_force.out');
            
            % Create .tcl file
            fid = fopen(filename_input,'w');
            fprintf(fid,'model BasicBuilder -ndm 1 -ndf 1 \n');
            for i = 0:obj.nStories
                fprintf(fid,'node %i 0.0\n',i);
            end
            for i = 1:obj.nStories
                fprintf(fid,'mass %i %g\n',i,obj.storyMass(i));
            end
            fprintf(fid,'fix 0 1 \n');
            for i = 1:length(obj.storySpringDefinition)
                fprintf(fid,'%s\n',obj.storySpringDefinition{i});
            end
            for i = 1:obj.nStories
                fprintf(fid,'element zeroLength %i %i %i -mat %i -dir 1\n',i,i-1,i,i);
            end
            
            fprintf(fid,'timeSeries Linear 1\n');
            fprintf(fid,'pattern Plain 1 1 {\n');
            for i = 1:obj.nStories
                if F(i) ~= 0
                    fprintf(fid,'    load %i %g \n',i,F(i));
                end
            end
            fprintf(fid,'} \n');
            
            fprintf(fid,'recorder Node -file {%s} -time -nodeRange 1 %i -dof 1 disp \n',filename_output_def,obj.nStories);
            fprintf(fid,'recorder Element -file {%s} -eleRange 1 %i force \n',filename_output_force,obj.nStories);
            fprintf(fid,'record \n');
            fprintf(fid,'system UmfPack \n');
            fprintf(fid,'constraints Penalty 1.0e12 1.0e12 \n');
            fprintf(fid,'test NormDispIncr 1e-5 10 1 \n');
            fprintf(fid,'algorithm Newton \n');
            fprintf(fid,'numberer RCM \n');
            fprintf(fid,'integrator DisplacementControl %i 1 %g\n',obj.controlStory1,obj.pushover_stepSize);
            fprintf(fid,'analysis Static \n');
            
            switch lower(type)
                case 'targetdrift'
                    fprintf(fid,'set ok [analyze %i]\n',ceil(targetDrift/obj.pushover_stepSize));
                    fprintf(fid,'if { $ok != 0 } {\n');
                    fprintf(fid,'    exit 1\n');
                    fprintf(fid,'}\n');
                case 'targetpostpeakratio'
                    fprintf(fid,'set currentLoad [getTime]\n');
                    fprintf(fid,'set maxLoad $currentLoad\n');
                    fprintf(fid,'while { $currentLoad >= [expr %g*$maxLoad] } {\n',targetPostPeakRatio);
                    fprintf(fid,'    set ok [analyze 1]\n');
                    fprintf(fid,'    if { $ok != 0 } {\n');
                    fprintf(fid,'        exit 2\n');
                    fprintf(fid,'    }\n');
                    fprintf(fid,'    set currentLoad [getTime]\n');
                    fprintf(fid,'    if { $currentLoad > $maxLoad } {\n');
                    fprintf(fid,'        set maxLoad $currentLoad\n');
                    fprintf(fid,'    }\n');
                    fprintf(fid,'    if { [nodeDisp %i 1] > %g } {\n',obj.controlStory1,obj.pushover_maxDrift);
                    fprintf(fid,'        exit 3\n');
                    fprintf(fid,'    }\n');                    
                    fprintf(fid,'}\n');                    
                otherwise
                    error('Unknown analysis type: %s',type);
            end

            fprintf(fid,'exit 1 \n');
            fclose(fid);
            
            % Run OpenSees
            [status, result] = obj.runOpenSees(filename_input);
            results.textOutput = result;
            switch status
                case 1
                    results.exitStatus = 'Analysis Successful';
                case 2
                    results.exitStatus = 'Analysis Failed';
                case 3
                    results.exitStatus = 'Peak Drift Reached';
                otherwise
                    fprintf('%s\n',result);
                    error('Analysis Failed in Unknown Manner (exit code: %i)',status);
            end
            
            % Read Results
            temp = dlmread(filename_output_def);
            time = temp(:,1);
            results.totalDrift = temp(:,2:end);
            temp = dlmread(filename_output_force);
            results.storyShear = temp(:,2:2:end);
            
            % Computed Results
            storyDrift = results.totalDrift;
            storyDrift(:,2:end) = storyDrift(:,2:end)-storyDrift(:,1:(end-1));
            results.storyDrift = storyDrift;            
            results.appliedStoryForce = time*F;
            results.roofDrift = results.totalDrift(:,end);
            results.baseShear = results.storyShear(:,1);
            
            % Clean Folder
            if obj.deleteFilesAfterAnalysis
                delete(filename_input,filename_output_def,filename_output_force);
            end
        end   
        function results = responseHistory(obj,groundMotionFilename,dt,SF,tend)
            
            % Initilize Results
            results = struct;
            
            % Filenames
            filename_input              = obj.scratchFile('mdofShearBuilding2d_input.tcl');
            filename_output_timeSeries  = obj.scratchFile('mdofShearBuilding2d_timeSeries.out');
            filename_output_def         = obj.scratchFile('mdofShearBuilding2d_disp.out');
            filename_output_force       = obj.scratchFile('mdofShearBuilding2d_force.out');
            
            % Create .tcl file
            fid = fopen(filename_input,'w');
            fprintf(fid,'model BasicBuilder -ndm 1 -ndf 1 \n');
            
            writeFunction_updateRayleighDamping(fid)
            
            for i = 0:obj.nStories
                fprintf(fid,'node %i 0.0\n',i);
            end
            for i = 1:obj.nStories
                fprintf(fid,'mass %i %g\n',i,obj.storyMass(i));
            end
            fprintf(fid,'fix 0 1 \n');
            for i = 1:length(obj.storySpringDefinition)
                fprintf(fid,'%s\n',obj.storySpringDefinition{i});
            end
            for i = 1:obj.nStories
                fprintf(fid,'element zeroLength %i %i %i -mat %i -dir 1\n',i,i-1,i,i);
            end
            
            
            fprintf(fid,'timeSeries Path 1 -dt %g -filePath {%s} -factor %g\n',dt,groundMotionFilename,SF);
            fprintf(fid,'pattern UniformExcitation 1 1 -accel 1\n');

            fprintf(fid,'recorder Node -file {%s} -timeSeries 1 -node 0 -dof 1 accel\n',filename_output_timeSeries);           
            fprintf(fid,'recorder Node -file {%s} -time -nodeRange 1 %i -dof 1 disp \n',filename_output_def,obj.nStories);
            fprintf(fid,'recorder Element -file {%s} -eleRange 1 %i force \n',filename_output_force,obj.nStories);
            fprintf(fid,'record \n');
            
            fprintf(fid,'system UmfPack \n');
            fprintf(fid,'constraints Transformation \n');
            fprintf(fid,'test NormDispIncr 1e-5 10 1 \n');
            fprintf(fid,'algorithm Newton \n');
            fprintf(fid,'numberer RCM \n');
            
            fprintf(fid,'updateRayleighDamping %i %g %i %g\n',...
                obj.damping_ModeA,obj.damping_RatioA,...
                obj.damping_ModeB,obj.damping_RatioB);
            
            fprintf(fid,'integrator Newmark 0.50 0.25\n');
            fprintf(fid,'analysis VariableTransient \n');
            
            fprintf(fid,'set currentTime [getTime]\n');
            fprintf(fid,'while { $currentTime < %g } {\n',tend);
            fprintf(fid,'    set ok [analyze 1 %g]\n',dt);
            fprintf(fid,'    if { $ok != 0 } {\n');
            fprintf(fid,'        exit 2\n');
            fprintf(fid,'    }\n');
            fprintf(fid,'    set currentTime [getTime]\n');                   
            fprintf(fid,'}\n');                    

            fprintf(fid,'exit 1 \n');
            fclose(fid);
            
            % Run OpenSees
            [status, result] = obj.runOpenSees(filename_input);
            results.textOutput = result;
            switch status
                case 1
                    results.exitStatus = 'Analysis Successful';
                case 2
                    results.exitStatus = 'Analysis Failed';
                otherwise
                    fprintf('%s\n',result);
                    error('Analysis Failed in Unknown Manner (exit code: %i)',status);
            end
            
            % Read Results
            temp = dlmread(filename_output_timeSeries);
            results.groundMotion = temp;
            temp = dlmread(filename_output_def);
            results.time = temp(:,1);
            results.totalDrift = temp(:,2:end);
            temp = dlmread(filename_output_force);
            results.storyShear = temp(:,2:2:end);
            
            % Computed Results
            storyDrift = results.totalDrift;
            storyDrift(:,2:end) = storyDrift(:,2:end)-storyDrift(:,1:(end-1));
            results.storyDrift = storyDrift;            
            results.roofDrift = results.totalDrift(:,end);
            results.baseShear = results.storyShear(:,1);
            
            % Clean Folder
            if obj.deleteFilesAfterAnalysis
                delete(filename_input,filename_output_timeSeries,...
                    filename_output_def,filename_output_force);
            end
        end        
    end
end


function tf = isvectorsize(v,n)
s = size(v);
tf = isequal(s,[n 1]) || isequal(s,[1 n]);
end

function writeFunction_updateRayleighDamping(fid)
fprintf(fid,'proc updateRayleighDamping { modeA ratioA modeB ratioB } {\n');
fprintf(fid,'    # ###################################################################\n');
fprintf(fid,'    # updateRayleighDamping $modeA $ratioA $modeB $ratioB\n');
fprintf(fid,'    # ###################################################################\n');
fprintf(fid,'    # Runs an eigenvalue analysis and set proportional damping based on\n');
fprintf(fid,'    # the current state of the structure\n');
fprintf(fid,'    #\n');
fprintf(fid,'    # Input Parameters:\n');
fprintf(fid,'    # modeA, modeB - modes that will have perscribed damping ratios\n');
fprintf(fid,'    # ratioA, ratioB - damping ratios perscribed at the specified modes\n');
fprintf(fid,'\n');
fprintf(fid,'    # Get natural frequencies at the desired modes\n');
fprintf(fid,'    if { $modeA > $modeB } {\n');
fprintf(fid,'        set maxMode $modeA\n');
fprintf(fid,'    } else {\n');
fprintf(fid,'        set maxMode $modeB\n');
fprintf(fid,'    }\n');
fprintf(fid,'\n');
fprintf(fid,'    set eigs    [eigen $maxMode]\n');
fprintf(fid,'    set freqA   [expr sqrt([lindex $eigs [expr $modeA-1]])]\n');
fprintf(fid,'    set freqB   [expr sqrt([lindex $eigs [expr $modeB-1]])]\n');
fprintf(fid,'\n');
fprintf(fid,'    # Compute the damping factors\n');
fprintf(fid,'    set tempVal [expr 2.0/($freqA*$freqA-$freqB*$freqB)]\n');
fprintf(fid,'    set aM      [expr $tempVal*$freqA*$freqB*($ratioB*$freqA-$ratioA*$freqB)]\n');
fprintf(fid,'    set aK      [expr $tempVal*($ratioA*$freqA-$ratioB*$freqB)]\n');
fprintf(fid,'\n');
fprintf(fid,'    # Set the damping\n');
fprintf(fid,'    rayleigh $aM 0.0 0.0 $aK\n');
fprintf(fid,'}\n');
end