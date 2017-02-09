classdef elasticAnalysis2d_equalDof
    %elasticAnalysis2d_equalDof Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        retainedJoint     % ID of Retained Joint
        constrainedJoint  % ID of Constrained Joint
        dofs = [0 0 0];   % Flag if dof is fixed or free (0 = free, 1 = fixed)
    end
    
    methods
        %% Constructor 
        function obj=elasticAnalysis2d_equalDof(retainedJoint,constrainedJoint,dofs)
            % elasticAnalysis2d_equalDof(retainedJoint,constrainedJoint,dofs)
            obj.retainedJoint = retainedJoint;
            obj.constrainedJoint = constrainedJoint;
            obj.dofs = dofs;
        end
        
        %% Set Functions
        function obj = set.retainedJoint(obj,retainedJoint)
            if isscalar(retainedJoint) && isnumeric(retainedJoint)
                obj.retainedJoint = retainedJoint;
            else
                error('retainedJoint should be a numeric scalar');
            end
        end
        function obj = set.constrainedJoint(obj,constrainedJoint)
            if isscalar(constrainedJoint) && isnumeric(constrainedJoint)
                obj.constrainedJoint = constrainedJoint;
            else
                error('constrainedJoint should be a numeric scalar');
            end
        end
        function obj = set.dofs(obj,dofs)
            % @todo - check that each element is either a zero or a one
            if isequal(size(dofs),[3 1])
                obj.dofs = dofs';
            elseif isequal(size(dofs),[1 3])
                obj.dofs = dofs;
            else
                error('dofs is the wrong size')
            end            
        end
        
        %% Output Functions
        function print(obj,fid,flag)
            % print()
            % print(fid)
            % print(fid,flag)
            if (nargin < 2); fid = 1; end
            if (nargin < 3); flag = 1; end

            switch flag
                case 1
                    fprintf(fid,'    Equal DOF Contraint: Retained Joint: %i, Retained Joint: %i, DOFs: %i %i %i\n',...
                        obj.retainedJoint,obj.constrainedJoint,obj.dofs);
                otherwise 
                    warning('elasticAnalysis2d:badInput','Unknown print flag')
            end
        end 
    end
    
end

