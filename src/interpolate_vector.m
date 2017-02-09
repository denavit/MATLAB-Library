function val = interpolate_vector(vec,ind,x)

assert(isvector(vec) && isnumeric(vec),'vec should be a numeric vector');
assert(~isempty(ind) && ~isempty(x),'invalid ind or x (empty)');
assert(ind > 0 && ind <= numel(vec),'invalid ind');

if ind < numel(vec)
    val = vec(ind) + x*(vec(ind+1)-vec(ind));
else
    assert(x == 0,'if ind < numel(vec) then x must be zero')
    val = vec(ind);
end
    
end