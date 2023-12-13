classdef interactiveMapVisualization < handle
    % DoA Visualization of Impulse Responses
    % Created by: Nils Meyer-Kahlen
    % Adapted by Aaron Geldert
    % Last modified: 10 Oct 2022
    
    properties
        % data
        irOmni
        doa
        
        fs
        irLen_samp
        tax
        
        % parameters
        irDyn
        mapDyn
        maxDotRadius
        printFileName
        plotname
        
        % current selection
        currectTimeSelection_ms
        currentDynLimit_db
        selected_samp
        
        name
        
        doaTrue
        radiusTrue
    end
    
    properties (Access = private)
        % handles
        fig
        s1
        s2
        hbox
        xCircle
        yCircle
        paCircle
        savebtn
        
    end
    
    methods
        function obj = interactiveMapVisualization(p, doa, pars)
            % INPUT: p ... len x 1 omni input response
            %        doa ... the direction of arrival data
            %        pars
            
            assert(size(doa, 1) == size(p, 1), ...
                'len x 3 doa data and omni ir must be the same length');
            
            % Load the IR
            obj.irOmni = p;
            
            % Load the doa
            obj.doa = doa;
                   
            obj.fs = pars.fs;
            
            obj.mapDyn = pars.dBDynamics;
            obj.irDyn = pars.dBDynamics;
            
            obj.name = pars.name;
            obj.printFileName = pars.name;
            
            normalize = false;
            
            if normalize
                obj.irOmni =  obj.irOmni / max(abs( obj.irOmni (:, 1)));
            end
            
            obj.irLen_samp = size(obj.irOmni, 1);
            obj.tax = (0:obj.irLen_samp-1) / obj.fs * 1000;

            % setup the figure itself
            obj.plotname = obj.name;
            obj.setupFigure();
            
            % initial plot of the IR
            obj.updateIrPlot();
            
            % setup the rectangle
            obj.hbox = drawrectangle(obj.s2);
            addlistener(obj.hbox,'ROIMoved',@obj.boxMoved);
            
            phiCircle=linspace(0,2*pi);
            obj.xCircle=cos(phiCircle);
            obj.yCircle=sin(phiCircle);
            
            % this can be used to add infomation to the plot,
            % as for example expected directions of reflections
            if isfield(pars, 'doaTrue')
                obj.doaTrue = pars.doaTrue;
                obj.radiusTrue = pars.radiusTrue;
            end
            
        end
        
        function [] = setupFigure(obj)
            obj.fig = figure;
            obj.s1 = subplot(3, 1, 1:2);
            obj.drawGrid();
            axis off
            axis([-1 1 -0.5 .5])
            obj.s2 = subplot(3, 1, 3);
            
            obj.savebtn = uicontrol('Style', 'pushbutton', 'String', 'Save and Close',...
                'Position', [0 0 100 30 ],...
                'Callback', @obj.savePlot);
            
            title(obj.s1, obj.plotname)
            
        end
        
        function [] = boxMoved(obj, ~, evt)
            boxPosition = evt.CurrentPosition;
            obj.currectTimeSelection_ms(1) = boxPosition(1);
            obj.currectTimeSelection_ms(2) = boxPosition(1)+boxPosition(3);
            
            obj.currentDynLimit_db =  boxPosition(2);
            
            %disp(['Dymamics ' num2str(obj.currentDynLimit_db) ' dB']);
            
            obj.updateMap();
            
        end
        
        function [] = updateMap(obj)
            
            obj.selected_samp = ceil(obj.currectTimeSelection_ms / 1000 * obj.fs);
            
            obj.selected_samp(1) = max(obj.selected_samp(1), 1);
            obj.selected_samp(2) = min(obj.selected_samp(2), obj.irLen_samp);
            
            ir_dB = obj.lin2db(obj.irOmni, obj.irDyn);
            
            delete(obj.paCircle);
            color = [0, 0, 0];
            
            % figure(100), plot(ir_dB), hold on
            % plot(obj.selected_samp(1):obj.selected_samp(2), ir_dB(obj.selected_samp(1):obj.selected_samp(2)))
            % hold off
            
            mapDynLimit = obj.irDyn - obj.mapDyn;
            
            for iDraw = obj.selected_samp(1):obj.selected_samp(2)
                if (ir_dB(iDraw)>=obj.currentDynLimit_db && ir_dB(iDraw)>=mapDynLimit)
                    
                    aziEleRad = cart2sphExpand(obj.doa(iDraw, :));
                    
                    xy = obj.project(aziEleRad(:, 1:2));
                    r = (ir_dB(iDraw) / obj.irDyn).^1.5 * 0.035;
                    
                    obj.paCircle = [obj.paCircle, ...
                        patch(obj.s1, xy(:, 1)+r*obj.xCircle,xy(:, 2)+r*obj.yCircle, 1, 'faceColor', color,...
                        'facealpha', 0.8, 'EdgeColor', 'none')];
                end
            end
            
            for iDoaTrue = 1:size(obj.doaTrue, 1)
                   aziEleRad = cart2sphExpand(obj.doaTrue(iDoaTrue, :));
                    
                    xy = obj.project(aziEleRad(:, 1:2));
                    r = obj.radiusTrue(iDoaTrue);
                    
                    obj.paCircle = [obj.paCircle, ...
                        patch(obj.s1, xy(:, 1)+r*obj.xCircle,xy(:, 2)+r*obj.yCircle, 1, 'faceColor', [1 0 0],...
                        'facealpha', 0.8, 'EdgeColor', 'none')];
            end
            
        end
        
        function [] = updateIrPlot(obj)
            
            ir_db = obj.lin2db(obj.irOmni(:, 1), obj.irDyn);
            plot(obj.s2, obj.tax, ir_db, 'k')
            obj.s2.YTickLabel =  - max(obj.s2.YTick) + obj.s2.YTick;
            
            xlabel(obj.s2, 'Time in ms')
            ylabel(obj.s2,'Amplitude in dB')
            grid(obj.s2, 'on')
            
        end
        
        function [x_db] = lin2db(~, x, dyn)
            
            x_db = max(db(x) + dyn, 0);
        end
        
        function [] = drawGrid(obj)
            
            aziGrid = [0 30  60 90 120 150 180 -30 -60  -120 -150 -90 179.9]'/180*pi;
            eleGrid = [-60 -30 0 30 60 ]'/180*pi;
            linegray = 0.7;
            
            for k=1:length(aziGrid)
                ele=linspace(-pi/2, pi/2, 50);
                azi=ones(size(ele))*aziGrid(k);
                xyGrid = obj.project([azi', ele']);
                plot(obj.s1, xyGrid(:, 1),xyGrid(:, 2),'--', 'color', [1 1 1]*linegray);
                hold on
                text(xyGrid(25, 1),0, num2str(round(180/pi*aziGrid(k))), ...
                    'fontsize', 8, 'color', [1 1 1]*linegray )
            end
            
            for k=1:length(eleGrid)
                azi=linspace(-pi, pi-.01, 50);
                ele=ones(size(azi))*eleGrid(k);
                xyGrid = obj.project([azi', ele']);
                plot(obj.s1, xyGrid(:, 1),xyGrid(:, 2),'--', 'color', [1 1 1]*linegray);
                hold on
                text(0, xyGrid(25, 2), num2str(180/pi*eleGrid(k)), ...
                    'fontsize', 8, 'color', [1 1 1]*linegray )
            end
        end
        
        function [] = savePlot(obj, src, ~, ~)
            
            ir_db = obj.lin2db(obj.irOmni, obj.irDyn);
            
            irMarked = ir_db(obj.selected_samp(1):obj.selected_samp(2));
            irMarked(irMarked <= obj.currentDynLimit_db) = nan;
            
            plot(obj.s2, obj.tax, ir_db, 'color', [1 1 1]*.8)
            hold on
            
            hstem = stem(obj.s2, obj.tax(obj.selected_samp(1):obj.selected_samp(2)), ...
                irMarked);
            
            hstem.Marker = 'none';
            hstem.Color = 'k'
            
            xlabel('Time (ms)')
            ylabel('Amplitude (dB)')
            grid('on')
            
            delete(src)
            delete(obj.hbox)
            obj.fig
            printScaled(14, 10, obj.printFileName, 'png')
            close(obj.fig)
           
        end
        
        function [xy] = project(obj, aziEle)
            % Hammer-Aidhof projection for now. Could be changed
            
            azi = aziEle(:, 1);
            ele = aziEle(:, 2);
            
            hap = @(azi, ele) [-cos(ele).*sin(azi/2) 0.5 * sin(ele)] ./ ...
                (sqrt(1+cos(ele).*cos(azi/2)));
            
            xy = hap (mod (azi+ pi, 2 * pi) - pi, ele);
            
        end
    end
end

function [aer] = cart2sphExpand(xyz)
% convert cartesian to spherical coordinates
% without having to index x y and z separetely
% nmk20

if size(xyz, 2) ~= 3
    xyz = xyz';
end

[azi ele r] = cart2sph(xyz(:, 1), xyz(:, 2), xyz(:, 3));
aer = [azi ele r];

end


