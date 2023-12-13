classdef IsmSimulator < handle
    %ISMSIMULATOR ImageSourceModelSimulator to generate image sources
    %   Based on Lauri Savioja's work (www.interactiveacoustics.com)
    % Created by: Aaron Geldert & Nils Meyer-Kahlen
    % Last modified: 16 Oct 2022
    
    properties
        fs = 48e3
        c = 343
        sourcePosition
        receiverPosition
        maxIsOrder = 3
        irLength_s = 0.2
        geomFile = 'cuboid.txt';
        
        % walls
        rCoeff = 0.7071
        
        % reciever
        
        % ISM simulator
        geometry
        furthestPt
        totalSurfaceArea
        imageSources
        visibleImageSources
        numVisibleIs = 0
        isPositions
        PC
        
        verbose = 1
        
    end
    
    methods
        function obj = IsmSimulator(geomFile)
            %ISMSIMULATOR Construct an instance of this class
            %   Detailed explanation goes here
            if nargin < 1
                geomFile = 'cuboid.txt';
                warning('No geometry file specified, loading default (cuboid.txt)');
            end
            obj.loadGeometry(geomFile);
        end
        
        function computeImageSources(obj, source, maxOrder)
            %COMPUTEIMAGESOURCES Recursive generation of IS tree
            %   refer to the interactiveAcoustics pseudocode
            
            numSurfaces = length(obj.geometry);
            for iSurface = 1:numSurfaces
                surface = obj.geometry(iSurface);
                % check that surface is not previous reflector
                if surface.index == source.previousReflector
                    continue;
                end
                
                % check that new source is further than parent from receiver
                
                % make sure that prev reflector and surface are facing
                
%                 if (source.previousReflector == -1) || (dot(obj.geometry(source.previousReflector).normal, obj.geometry(surface.index).normal) <= 0)
                    % make sure that the surface is somewhat in front
                    % of the previous reflector - ALWAYS TRUE W/CONVEX GEOM
%                     if (source.previousReflector == -1) || (obj.geometry(source.previousReflector).inFront(obj.geometry(surface.index)))
                        % now we reflect the source against surface!

                        newSource = reflectSource(source, surface);
                        try
                            if isnan(newSource)
                                continue;
                            end
                        catch
                        end
                        
                        newSource.reflector = surface;
                        newSource.previousReflector = surface.index;
                        newSource.parent = source;
                        newSource.order = source.order + 1;

                        newSource.reflectors = [newSource.reflectors; newSource.parent.reflectors; surface];
                        obj.imageSources = [obj.imageSources; newSource];
                        if maxOrder > 0
                            % recursion, limited by maxOrder
                            obj.computeImageSources(newSource, maxOrder-1)
                        end
%                     end
%                 end
            end
            
            
            
            end
            
            function path = constructImageSourcePath(obj, is_, listenerPosition)
                % ref: http://interactiveacoustics.info/html/GA_IS_isFull.html
                
                originalIs = is_;
                originalIs.valid = 1;
                
                path = NaN(is_.order+1,3);
                path(is_.order+1,:) = listenerPosition;
                
                for order = originalIs.order:-1:1
                    intersectionPoint = is_.reflector.doesLineIntersect(is_.position, path(order+1,:));
                    if isscalar(intersectionPoint) && intersectionPoint == 0
                        originalIs.valid = 0;
                        path = 0; return;
                    end
                    path(order,:) = intersectionPoint;
                    is_ = is_.parent;
                end
                
                path(1,:) = is_.position;
                originalIs.path = path;
                originalIs.valid = 1;
            end
            
            % IGNORED IN CURRENT IMPLEMENTATION!
            function validateImageSource(obj, is_)
                % check if IS is occluded by some obstructing surfaces
                
                 if is_.valid
                    for order = 1:(is_.order)
                        segmentStart = is_.path(order,:);
                        segmentEnd = is_.path(order+1,:);
                        
                        reflector = is_.reflectors(order);
                        prevReflector = [];
                        if order > 1
                            prevReflector = is_.reflectors(order);
                        end
                        
                        for iSurface = 1:length(obj.geometry)
                            surface = obj.geometry(iSurface);
                            if surface.index ~= reflector.index
                            	if(~isempty(prevReflector) && (surface.index == prevReflector.index))
                                    continue;
                                end
%                                 disp(surface.index);
                                intersectionPoint = surface.doesLineIntersect(segmentStart, segmentEnd);
                                if ~isscalar(intersectionPoint)
                                    is_.valid = 0; return;
                                end
                            end
                        end
                    end
                end
            end
            
            function getVisibleSources(obj, listenerPosition)
                obj.visibleImageSources = [];
                numSources = length(obj.imageSources);
                for iSource = 1:numSources
                    source = obj.imageSources(iSource);
                    obj.constructImageSourcePath(source, listenerPosition);
                    % ignore occluded IS problems for now!
%                     obj.validateImageSource(source);
                    if source.valid
                        obj.visibleImageSources = [obj.visibleImageSources; source];
                    end
                end
            end
        
            function simulateImageSources(obj)
                % determines the image sources, but not visibility
                % follow up with imageSourcesForRcvPos to get PC of IS
                obj.imageSources = [];
                source = Source(obj.sourcePosition);
                obj.imageSources = source;
                obj.computeImageSources(source, obj.maxIsOrder-1);
                if obj.verbose ~= 0
                    disp(['Calculated ' num2str(length(obj.imageSources)) ' potential image sources.']);    
                end
            end
            
            function PC = imageSourcesForRcvPos(obj, rcvPos)
                % this part depends on rcv position
                obj.receiverPosition = rcvPos;
                obj.getVisibleSources(obj.receiverPosition);
                obj.numVisibleIs = length(obj.visibleImageSources);
                
                obj.isPositions = NaN(obj.numVisibleIs, 3);
                orders = NaN(obj.numVisibleIs,1);
                
                for ii = 1:obj.numVisibleIs
                    obj.isPositions(ii,:) = obj.visibleImageSources(ii).position - obj.receiverPosition;
                    orders(ii) = obj.visibleImageSources(ii).order;
                end
                
                PC.pos = obj.isPositions; % these are relative to receiver!
                PC.mass = ones(obj.numVisibleIs,1).*(obj.rCoeff.^orders); % -3 dB per reflection
                PC.mass = PC.mass ./ (vecnorm(PC.pos,2,2)+0.01); % 1/r attenuation with a regularization

                PC.n = obj.numVisibleIs;
                if obj.verbose ~= 0
                    disp(['Generated PC with ' num2str(PC.n) ' image sources visible.']);
                end
                obj.PC = PC;
            end
            
            function [hs, hr, hi] = plotGeometry(obj, showVS, showNormals)
                if nargin < 3
                    showNormals = 0;
                    if nargin < 2
                        showVS = 1;
                    end
                end
                hs = NaN; % initialize
                hr = NaN; 
                hi = NaN;
                
                strength = 20;
                if ~isempty(obj.receiverPosition)
                    if showVS ~= 0
                        hi = scatter3(obj.PC.pos(:,1) + obj.receiverPosition(1),...
                            obj.PC.pos(:,2) + obj.receiverPosition(2),...
                            obj.PC.pos(:,3) + obj.receiverPosition(3), 1+db(1+strength*obj.PC.mass).^1.4, 'o', 'filled', 'MarkerFaceColor', '#7E2F8E');
                        hold on;
                    end
                    hr = scatter3(obj.receiverPosition(1), obj.receiverPosition(2), obj.receiverPosition(3),...
                        50, 'g^', 'filled', 'MarkerEdgeColor', 'k');
                    hold on;
                    % mark the receiver relative to the ground
                    line([obj.receiverPosition(1) obj.receiverPosition(1)],...
                        [obj.receiverPosition(2) obj.receiverPosition(2)],...
                        [0 obj.receiverPosition(3)],'Color','black','LineStyle',':');
                    scatter3(obj.receiverPosition(1), obj.receiverPosition(2), 0,...
                        'k.');
                    hold on;
                end
                
                numSurfaces = length(obj.geometry);
                for iSurface = 1:numSurfaces
                    surface = obj.geometry(iSurface);
                    pts = surface.points;
                    
                    % horizontal planes plotted with no Edge Alpha
                    if length(unique(surface.points(:,3)))<=1
                        patch(pts(:,1), pts(:,2), pts(:,3), [.5 .5 .5],'FaceAlpha',.2,'EdgeAlpha',0);
                    else
                        patch(pts(:,1), pts(:,2), pts(:,3), [.5 .5 .5],'FaceAlpha',.2);
                    end
                    hold on;
                    
                    % plot the normal vector
                    if showNormals ~= 0
                        hq = quiver3(surface.center(1), surface.center(2), surface.center(3),...
                            surface.normal(1), surface.normal(2), surface.normal(3));
                        hq.Marker = '.';
                        hq.AutoScaleFactor = 0.5;
                        hq.MaxHeadSize = 2;
                    end
                end
                
                if (~showVS && ~isempty(obj.sourcePosition))
                    hs = scatter3(obj.sourcePosition(1), obj.sourcePosition(2), obj.sourcePosition(3),...
                        50, 'rs', 'filled', 'MarkerEdgeColor', 'k');
                    hold on;
                    line([obj.sourcePosition(1) obj.sourcePosition(1)],...
                        [obj.sourcePosition(2) obj.sourcePosition(2)],...
                        [0 obj.sourcePosition(3)],'Color','black','LineStyle',':');
                    scatter3(obj.sourcePosition(1), obj.sourcePosition(2), 0,...
                        'k.');
                end

                xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
                view(3);
                axis image;
                grid off;
                
            end
            
            function ret = isPointInsideGeometry(obj, pointInQuestion, offset)
                % assuming all surface normals point INSIDE the geometry
                % consider only points that are offset away from the wall
                if nargin < 3
                    offset = 0;
                end
                
                if any(isnan(pointInQuestion))
                    ret = 0; return;
                end
                for iSurface = 1:length(obj.geometry)
                    surface = obj.geometry(iSurface);
                    if dot(surface.normal, pointInQuestion-(surface.center + offset*surface.normal)) < 0
                        ret = 0; return;
                    end
                end
                ret = 1;
            end
            
            function loadGeometry(obj, filename)
                obj.geomFile = filename;
                iPoint = 0;
                numSurfaces = 0;
                pointsForSurface = zeros(4,3);
                surfaces = [];
                rectDims = zeros(1,3);
                obj.totalSurfaceArea = 0;

                lines = readlines(filename);
                for iLine = 1:length(lines)
                    % try to read coordinates, NaN if invalid
                    coords = str2double(split(lines(iLine),' '));
                    if any(isnan(coords))
                        continue;
                    else
                        % update furthest point (yes, per x y z dimension!)
                        rectDims = max(coords.', rectDims);

                        iPoint = iPoint + 1;
                        pointsForSurface(iPoint,:) = coords.';
                        if iPoint == 4
                            % add a surface
                            newSurface = Surface(pointsForSurface);
                            surfaces = [surfaces; newSurface];
                            numSurfaces = numSurfaces + 1;
                            iPoint = 0;
                        end
                    end
                end
                if isempty(surfaces)
                    warning(['No geometry data could be loaded from ' filename]);
                else
                    disp(['Loaded ' num2str(numSurfaces) ' surfaces from ' filename]);
                    for iSurface = 1:numSurfaces
                        obj.totalSurfaceArea = obj.totalSurfaceArea + surfaces(iSurface).area;
                    end
                    
                    obj.furthestPt = rectDims;
                    obj.geometry = surfaces;
                    
                end
                
            end
    end
end

