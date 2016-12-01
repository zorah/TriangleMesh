function [ dist, path ] = dijkstra_geodesic( M, a, b )
%DIJKSTRA_GEODESIC Calculates the distance between a and b or the distance map of a.
% Inputs:
%   a, b - indices of vertices of M
%   M - TriangleMesh
% Output:
%   dist - scalar if a and b were given, n x 1 if only a was given
%   path - 0 if the output is a distance map, otherwise array of indices
%
%   copyright (c) 2016 Zorah Lähner (laehner@in.tum.de)

if (~exist('b', 'var'))
    b = 0;
end

n = size(M.VERT,1);
m = size(M.TRIV,1);

% calculate length of all edges + adjacency list
mesh = cell(n,1); % mesh{i} is a 2xk matrix where k is number of neighbors of vertex i
mesh(:) = {[]};

for k=1:m
    indices = [M.TRIV(k,:), M.TRIV(k,1)];
    for i=1:3
        i1 = indices(i);
        i2 = indices(i+1);
        length = sqrt(sum((M.VERT(i1,:) - M.VERT(i2,:)).^2));
        mesh{i1} = [mesh{i1}, [i2; length]];
    end
end

% init heap
heap = MinHeap(n);

heap.push(0, a);

for i=1:n
    heap.push(Inf, i);
end

% run dijkstra
result = -1*ones(n,1);
predecessors = zeros(n,1);

while(~heap.isEmpty())
    [heap, current] = heap.pop();
    result(current(1,2)) = current(1,1);
    
    % update all neighbors of the current vertex
    neighbors = mesh{current(1,2)};
    for i=1:size(neighbors,2)
        index = neighbors(1,i);
        if(result(index) < 0) % path not minimal yet
            oldValue = heap.peakKey(index);
            newValue = current(1,1) + neighbors(2,i);
            
            % update?
            if(newValue < oldValue)
                % update!
                heap.decrease(index, newValue);
                predecessors(index) = current(1,2);
            end
        end
    end
    if (current(1,2) == b)
        break;
    end
end

if b ~= 0
    dist = result(b);
    path = b;
    pred = b;
    while(pred ~= a)
        pred = predecessors(pred);
        path = [path, pred];
        if pred < 1
            break;
        end
    end
    path = [path, a];
else
    dist = result;
    path = 0;
end

end

