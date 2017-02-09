function [ind,x] = find_limit_point_in_vector(vec,val)

assert(isvector(vec) && isnumeric(vec),'vec should be a numeric vector')
assert(isscalar(val) && isnumeric(val),'val should be a numeric scalar')

if vec(1) == val
    ind = 1;
    x = 0.0;
    return
elseif vec(1) > val
    vec2 = vec <= val;
else % vec(1) < val
    vec2 = vec >= val;
end

if isequal(zeros(size(vec2)),vec2)
    ind = [];
    x = [];
else
    ind = find(vec2, 1, 'first');
    if val == vec(ind)
        x = 0;
    else
        ind = ind-1;
        x = (val-vec(ind))/(vec(ind+1)-vec(ind));
    end
end
    
end