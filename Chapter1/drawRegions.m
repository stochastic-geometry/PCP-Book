% DRAWREGIONS   Voronoi graphical realization
%   This function will take in a set of generator points and regions and
%   will draw the associated graph
%   
%   ***** Inputs *****
%   bounds - the bounding positively oriented polygon
%   regionData - the regions are a nx3 cell array with the {i,1} element 
%       containing generator point information, the {i,2} element 
%       containing the arc information, and the {i,3}
%       containing edge segments (from intersections with the boundary)
%
%   Author: Andrew Kwok, University of California at San Diego
%   Date: 17 August 2009
function drawRegions(bounds, regionData,varargin)
%clf
n = size(regionData, 1);

%hold on
% Define agent label position
delta = max(max(bounds(1,:)) - min(bounds(1,:)), ...
    max(bounds(2,:)) - min(bounds(2,:))) / 60;
for i = 1:n
    pt = regionData{i, 1};
    arcs = regionData{i, 2};
    edges = regionData{i, 3};
    
    % Draw arcs and edges
    for j = 1:size(arcs, 1)
        drawArc(arcs(j, 1:2), arcs(j, 3), arcs(j, 4), arcs(j, 5),varargin{:});
    end
    for j = 1:size(edges, 1)
        line(edges(j, [1, 3]), edges(j, [2, 4]),varargin{:});
    end
end

% Plot the bounding polygon
% bounds(:, end+1) = bounds(:,1);
% plot(bounds(1,:), bounds(2,:), 'r', 'LineWidth', 2);
% axis equal
% axis square
% xlim([min(bounds(1,:)), max(bounds(1,:))]);
% ylim([min(bounds(2,:)), max(bounds(2,:))]);
% axis off
% hold off