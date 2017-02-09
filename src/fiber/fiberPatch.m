classdef fiberPatch
    methods (Abstract)
        [m,A,z] = fiberData2d(obj,axis,sf)
        [m,A,z,y] = fiberData3d(obj,sfz,sfy)
        [zmin,zmax,ymin,ymax] = bounds(obj)
        mat = matIDs(obj)
        plot(obj,lineWidth)
    end
end
