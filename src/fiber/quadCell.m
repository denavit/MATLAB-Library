function [A,z,y] = quadCell(vertexCoords)
assert(isequal([4 2],size(vertexCoords)),'vertexCoords must be a 4x2 matrix')
vertexCoords = vertcat(vertexCoords,vertexCoords(1,:));
A = 0;
z = 0;
y = 0;
for i = 1:4
    temp = (vertexCoords(i,1)*vertexCoords(i+1,2)-...
        vertexCoords(i+1,1)*vertexCoords(i,2));
    A = A + temp;
    z = z + (vertexCoords(i,1) + vertexCoords(i+1,1))*temp;
    y = y + (vertexCoords(i,2) + vertexCoords(i+1,2))*temp;
end
A = -A/2;
z = -z/(6*A);
y = -y/(6*A);
end