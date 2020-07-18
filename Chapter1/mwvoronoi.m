% MWVORONOI   energy weighted multiplicatively-weighted Voronoi 
%   partition generator
%   MWVORONOI will take in a bounding region along with a set of 
%   generator points and produce the ordinary Voronoi partition defined by
%       1/(E_i)^2 |q - p_i|^2 <= 1/(E_j)^2 |q - p_j|^2
%   for all j ~= i.
%
%   ***** Inputs *****
%   bounds - a 2 x n_b array specifying the positively oriented bounding
%       polygon in 2D
%   points - the 3 x n_p array specifying generator points.  The (3,i)
%       value describes the energy state E_i of a generator point.
%
%   ***** Outputs *****
%   regionData - a n x 3 cell with:
%       -- the {i,1} element containing generator point information, 
%       -- the {i,2} element containing the arc information. These arcs
%          are oriented and each row represents an arc:
%          [x, y, r, th1, th2]
%          (x, y) is the center of the arc, r is the radius of the arc,
%          th1 is the starting angle, and th2 is the ending angle.
%          If th1 < th2, then the region lies to the left (or inside)
%          the arc.  If th1 > th2, then the region lies to the right (or
%          outside) the arc.
%       -- the {i, 3} element containing edge segments (from intersections
%          with the boundary). Again, these are oriented, and each row
%          represents a segment:
%          [x1, y1, x2, y2]
%          And the region is always to the left of the vector from (x1, y1)
%          to (x2, y2).
%
%   Author: Andrew Kwok, University of California at San Diego
%   Date: 17 August 2009
function regionData = mwvoronoi(bounds, points)
n = size(points, 2);
regionData = cell(n, 3);

% Deal with numerical issues when Ei == Ej
for i = 1:n
    Ei = points(3, i);
    for j = i+1:n
        Ej = points(3, j);
        if (abs(Ei / Ej - 1) < 0.0001)
            points(3, j) = Ej * 0.9999;
        end
    end
end

for i = 1:n
%     arcs = [points(1, i), points(2, i), points(3, i), 0, 2*pi];
    arcs = zeros(0, 5);
    Pi = [points(1, i); points(2, i)];
    Ei = points(3, i);
    
    % Now iterate over all other points to construct the MW Voronoi Diagram
    for j = 1:n
        if(j == i)
            continue;
        end
        
        Pj = [points(1, j); points(2, j)];
        Ej = points(3, j);
        
        % define the center and radius of the Appolonius circle
        d = Pj - Pi;
        C = Pi + Ei^2 / (Ei^2 - Ej^2) * d;
        R = abs(Ei * Ej * norm(d) / (Ei^2 - Ej^2));
        
        % Now intersect this with all arcs so far
        [arcIndex, arcAngle, appoAngle] = findIntersects(C, R, arcs);
        intersects = [arcIndex, arcAngle, appoAngle];
        
        if (size(intersects, 1) > 0)
            % This returns a sorted list of valid arc segments
            [appoArcs, intersectStart, intersectEnd] = ...
                getArcs(C, R, Ei, Ej, intersects, arcs);

            % Now insert these arcs into the existing set of arcs
            % The new arcs themselves divide each existing arc into
            % sections that are still within the region and segments that
            % are not.  So crop out the ones that are not.
            if (Ei < Ej)
                appoCircle = [C', R, 0, 2*pi];
            else
                appoCircle = [C', R, 2*pi, 0];
            end
            
            newArcs = zeros(0, 5);
            % For each arc that isnt affected by the new Appolonius circle,
            % decide if it stays or goes
            rows = setdiff(1:size(arcs, 1)', ...
                [intersectStart(:, 1); intersectEnd(:, 1)]);
            testArcs = arcs(rows, :);
            for k = 1:size(testArcs, 1)
%                 testTh = (testArcs(k, 4) + testArcs(k, 5)) / 2;
%                 testPt = testArcs(k, 1:2)' + testArcs(k, 3) * ...
%                     [cos(testTh); sin(testTh)];
                if (isArcInside(testArcs(k, :), appoCircle))
                    newArcs = [newArcs; testArcs(k, :)];
                end
            end
            
            % Now examine the affected arcs and decide which parts to keep
            intersects = [intersectStart; intersectEnd];
            rows = unique(intersects(:, 1));
            for k = 1:length(rows)
                indices = find(intersects(:, 1) == rows(k));
                angles = intersects(indices, 2);
                testArcs = getArcSegments(arcs(rows(k), :), angles);
                for l = 1:size(testArcs, 1)
%                     testTh = (testArcs(l, 4) + testArcs(l, 5)) / 2;
%                     testPt = testArcs(l, 1:2)' + testArcs(l, 3) * ...
%                         [cos(testTh); sin(testTh)];
                    if (isArcInside(testArcs(l, :), appoCircle))
                        newArcs = [newArcs; testArcs(l, :)];
                    end
                end
            end
            
            % Finally add the new Appolonius arcs
            for k = 1:length(appoArcs)
                newArcs = [newArcs; appoArcs{k}];
            end
        else % No intersection between existing region and new appo Circle
            if (size(arcs, 1) > 0)
                if (Ei < Ej)
                    appoCircle = [C', R, 0, 2*pi];
                else
                    appoCircle = [C', R, 2*pi, 0];
                end
                
                newArcs = zeros(0, 5);
                % Now figure out which arcs to keep and which to throw
                % away
                for k = 1:size(arcs, 1)
                    if (isArcInside(arcs(k, :), appoCircle))
                        newArcs = [newArcs; arcs(k, :)];
                    end
                end
                
                % And see if the new appo circle is to be included
                if (isArcInside(appoCircle, arcs))
                    newArcs = [newArcs; appoCircle];
                    
                    
                end
                
            else
                if (Ei < Ej)
                    newArcs = [C', R, 0, 2*pi];
                else
                    newArcs = [C', R, 2*pi, 0];
                end
            end
        end
        arcs = newArcs;
    end
    regionData{i, 1} = points(:, i);
    regionData{i, 2} = arcs;
end

% and now intersect what we have with the bounding polygon, Q
extendedBounds = [bounds, bounds(:, 1)];
for i = 1:n
    region = regionData{i, 2};
    newArcs = zeros(0, 5);
    edges = zeros(0, 4);
    
    % find the intersections of each boundary segment with one region
    % intersects has columns corresponding to
    % [regionArcIndex, segmentParameter, arcParameter]
    % segmentParameter is a number in [0, 1]
    % arcParameter is a number in [0, 2*pi]
    intersects = zeros(0, 3);
    for j = 1:size(bounds, 2)
        p1 = extendedBounds(:, j);
        p2 = extendedBounds(:, j+1);
        
        thisIntersects = findLineArcIntersects(p1, p2, region);
        % Determine which segment to include in the region (if any)
        if (size(thisIntersects, 1) > 0)
            segments = getSegments(p1, p2, thisIntersects(:, 2));
            for k = 1:size(segments, 1)
                testPt = (segments(k, 1:2)' + segments(k, 3:4)') / 2;
                if (isPtInside(testPt, region))
                    edges = [edges; segments(k, :)];
                end
            end
        else % check if the entire segment is good
            testPt = (p1 + p2) / 2;
            if (isPtInside(testPt, region))
                edges = [edges; [p1', p2']];
            end
        end
        
        intersects = [intersects; thisIntersects];
    end
    
    % For each arc in region that does not intersect bounds, decide if
    % it stays or goes
    rows = setdiff(1:size(region, 1), intersects(:, 1));
    testArcs = region(rows, :);
    for j = 1:size(testArcs, 1)
        testTh = (testArcs(j, 4) + testArcs(j, 5)) / 2;
        testPt = testArcs(j, 1:2)' + testArcs(j, 3) * ...
            [cos(testTh); sin(testTh)];
        if (isInBounds(testPt, bounds))
            newArcs = [newArcs; testArcs(j, :)];
        end
    end
    
    % Now examine the affected arcs and see which portions to keep
    rows = unique(intersects(:, 1));
    for j = 1:length(rows)
        indices = find(intersects(:, 1) == rows(j));
        angles = intersects(indices, 3);
        testArcs = getArcSegments(region(rows(j), :), angles);
        for k = 1:size(testArcs, 1)
            testTh = (testArcs(k, 4) + testArcs(k, 5)) / 2;
            testPt = testArcs(k, 1:2)' + testArcs(k, 3) * ...
                [cos(testTh); sin(testTh)];
            if (isInBounds(testPt, bounds))
                newArcs = [newArcs; testArcs(k, :)];
            end
        end
    end
    
    regionData{i, 2} = newArcs;
    regionData{i, 3} = edges;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Helper function to find intersection of lines with circles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function intersects = findLineArcIntersects(p1, p2, region)
d = p2 - p1;
% Iterate through each arc of the region to see if there is an intersection
intersects = zeros(0, 3);
for i = 1:size(region, 1)
    center = region(i, 1:2)';
    R = region(i, 3);
    a = -d(1)^2 - d(2)^2;
    b = 2 * (center(1)*d(1) + center(2)*d(2) - d(1)*p1(1) - d(2)*p1(2));
    c = -center(1)^2 - center(2)^2 + 2*center(1)*p1(1) - p1(1)^2 + ...
        2*center(2)*p1(2) - p1(2)^2 + R^2;
    if (b^2 - 4*a*c > 0) % Intersection exists
        tPlus = (-b + sqrt(b^2 - 4*a*c)) / (2*a);
        tMinus = (-b - sqrt(b^2 - 4*a*c)) / (2*a);
        if (0 < tMinus && tMinus < 1)
            vec = p1 + d*tMinus - center;
            theta = atan2(vec(2), vec(1));
            if (theta < 0)
                theta = theta + 2*pi;
            end
            alpha = min(region(i, 4), region(i, 5));
            beta = max(region(i, 4), region(i, 5));
            if (alpha <= theta && theta < beta)
                intersects(end+1, :) = [i, tMinus, theta];
            end
        end
        
        if (0 < tPlus && tPlus < 1)
            vec = p1 + d*tPlus - center;
            theta = atan2(vec(2), vec(1));
            if (theta < 0)
                theta = theta + 2*pi;
            end
            alpha = min(region(i, 4), region(i, 5));
            beta = max(region(i, 4), region(i, 5));
            if (alpha <= theta && theta < beta)
                intersects(end+1, :) = [i, tPlus, theta];
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Helper function to determine if a point lies inside a region
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function bool = isPtInside(point, region)
bool = 1;
for i = 1:size(region, 1)
    center = region(i, 1:2)';
    R = region(i, 3);
    alpha = region(i, 4);
    beta = region(i, 5);
    if (alpha < beta) % positively oriented, check if p1, p2 is inside
        if (norm(point - center) > R)
            bool = 0;
            return;
        end
    else % negatively oriented, check if p1, p2 is outside
        if (norm(point - center) < R)
            bool = 0;
            return;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Helper function to determine if a point lies inside a region
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function bool = isArcInside(arc, region)
bool = 1;

C = arc(1:2)';
R = arc(3);
alpha = arc(4);
beta = arc(5);

% Test the endpoints and the midpoint of the arc
testTh = alpha + 0.001 * (beta - alpha);
testPt = C + R * [cos(testTh); sin(testTh)];
if (~isPtInside(testPt, region))
    bool = 0;
    return;
end

testTh = alpha + 0.5 * (beta - alpha);
testPt = C + R * [cos(testTh); sin(testTh)];
if (~isPtInside(testPt, region))
    bool = 0;
    return;
end

testTh = alpha + 0.999 * (beta - alpha);
testPt = C + R * [cos(testTh); sin(testTh)];
if (~isPtInside(testPt, region))
    bool = 0;
    return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Helper function to determine if a point lies inside bounds
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function bool = isInBounds(point, bounds)
% Assuming bounds is positively oriented!
bounds = [bounds, bounds(:, 1)];
point = point(1:2);

bool = 1;
for i = 1:size(bounds, 2)-1
    boundVec = bounds(:, i+1) - bounds(:, i);
    ptVec = point - bounds(:, i);
    cross = boundVec(1) * ptVec(2) - ptVec(1) * boundVec(2);
    if (cross < 0)
        bool = 0;
        return;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Helper function to find the intersection of a circle with a set of arcs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [arcIndex, arcAngle, appoAngle] = findIntersects(C, R, arcs)
arcIndex = [];
arcAngle = [];
appoAngle = [];

for i = 1:size(arcs, 1)
    d = C - [arcs(i, 1); arcs(i, 2)];
    gamma = acos((arcs(i, 3)^2 - R^2 + norm(d)^2) / (2 * arcs(i, 3) * norm(d)));
    
    if (imag(gamma) == 0)   % There is an intersection
        th = atan2(d(2), d(1));
        alpha = th - gamma;
        beta = th + gamma;
        % alpha and beta need to be in [0, 2 pi]
        alpha = mod(alpha, 2*pi);
        beta = mod(beta, 2*pi);
        
        % Test to see if alpha or beta lie between the angles in arcs(i, :)
        arcAlpha = min(arcs(i, 4), arcs(i, 5));
        arcBeta = max(arcs(i, 4), arcs(i, 5));
        phi = acos((R^2 - arcs(i, 3)^2 + norm(d)^2) / (2 * R * norm(d)));
        if (arcAlpha < alpha && alpha < arcBeta)
            nextAppoAngle = mod(pi + th + phi, 2*pi);
            arcIndex = [arcIndex; i];
            arcAngle = [arcAngle; alpha];
            appoAngle = [appoAngle; nextAppoAngle];
        end
        if (arcAlpha < beta && beta < arcBeta)
            nextAppoAngle = mod(pi + th - phi, 2*pi);
            arcIndex = [arcIndex; i];
            arcAngle = [arcAngle; beta];
            appoAngle = [appoAngle; nextAppoAngle];
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Helper function to return the valid Appolonius arcs from the intersection
% of a Appolonius circle with an existing region
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [appoArcs, intersectStart, intersectEnd] = ...
    getArcs(C, R, Ei, Ej, intersects, arcs);
% Create a list of candidate arcs. If Ei > Ej, then the additional
% arcs must be negatively oriented
intersects = sortrows(intersects, 3);
if (Ei > Ej)
    intersects = flipud(intersects);
end
appoArcs = cell(size(intersects, 1), 1);
for k = 1:size(intersects, 1)-1
    appoArcs{k} = [C', R, intersects(k, 3), ...
        intersects(k+1, 3)];
end
if (Ei > Ej)
    appoArcs{end} = [C', R, intersects(end, 3), 0; ...
        C', R, 2*pi, intersects(1, 3)];
else
    appoArcs{end} = [C', R, intersects(end, 3), 2*pi; ...
        C', R, 0, intersects(1, 3)];
end


% Now test each arc to see if it lies inside the existing region.
realAppoArcs = cell(0, 1);
intersectStart = zeros(0, 3);
intersectEnd = zeros(0, 3);
for k = 1:size(appoArcs, 1);
    % testTh = (appoArcs{k}(1, 4) + appoArcs{k}(1, 5)) / 2;
%     testTh = appoArcs{k}(1, 4) + 0.001 * (appoArcs{k}(1, 5) - appoArcs{k}(1, 4));
%     testPt = C + R * [cos(testTh); sin(testTh)];
    if (isArcInside(appoArcs{k}(1,:), arcs))
        realAppoArcs{end+1} = appoArcs{k};
        intersectStart(end+1, :) = intersects(k, :);
        if (k == size(appoArcs, 1))
            intersectEnd(end+1, :) = intersects(1, :);
        else
            intersectEnd(end+1, :) = intersects(k+1, :);
        end
    end
end

% Now sort the appoArcs according to positive orientation
% along the existing region.
[appoArcs, intersectStart, intersectEnd] = ...
    sortArcs(arcs, realAppoArcs, intersectStart, intersectEnd);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Helper function to sort Appolonius arcs of intersection along an existing
% positively oriented region
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [appoArcs, intersectStart, intersectEnd] = ...
    sortArcs(arcs, appoArcs, intersectStart, intersectEnd)

intersects = [intersectStart, intersectEnd];
% The following sorts the appoArcs according to the arcIndex of arcs
[intersects, permutation] = sortrows(intersects, 1);
appoArcs = appoArcs(permutation)';
intersectStart = intersects(:, 1:3);
intersectEnd = intersects(:, 4:6);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Helper function to return the arc segments that result from the
% intersection of an arc at certain angles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function arcSegments = getArcSegments(arc, angles)
arcSegments = zeros(length(angles)+1, 5);
if (arc(4) < arc(5))
    angles = sort(angles, 'ascend');
else
    angles = sort(angles, 'descend');
end

arcSegments(1, 1:4) = arc(1:4);
arcSegments(1, 5) = angles(1);
for i = 2:size(arcSegments, 1)-1
    arcSegments(i, 1:3) = arc(1:3);
    arcSegments(i, 4) = angles(i-1);
    arcSegments(i, 5) = angles(i);
end
arcSegments(end, 1:3) = arc(1:3);
arcSegments(end, 4) = angles(end);
arcSegments(end, 5) = arc(5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Helper function to return segments that result from intersecting a line
% segment with a region
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function segments = getSegments(p1, p2, segParams)
segments = zeros(length(segParams)+1, 4);
segParams = sort(segParams, 'ascend');

startPt = p1;
endPt = p1 + segParams(1) * (p2 - p1);
segments(1, :) = [startPt', endPt'];
for i = 2:length(segParams)
    startPt = p1 + segParams(i-1) * (p2 - p1);
    endPt = p1 + segParams(i) * (p2 - p1);
    segments(i, :) = [startPt', endPt'];
end
startPt = p1 + segParams(end) * (p2 - p1);
endPt = p2;
segments(end, :) = [startPt', endPt'];