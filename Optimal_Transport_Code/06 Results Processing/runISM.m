function [PCs, config] = runISM(params)
% Created by: Aaron Geldert
% Last modified: 10 Oct 2022
    
    switch params.moving
        case 'rcv'
            config.movDir = params.movDir;
            config.movInterval = params.movInterval;
            config.numMovPos = params.numMovPos;
            movStep = params.movInterval .* params.movDir./vecnorm(params.movDir,2);
            
            config.srcPos = params.srcPos;
            config.rcvPos = ones(params.numMovPos,3).*params.rcvPos +...
                (0:(params.numMovPos-1)).' * movStep;

        case 'src'
            config.movDir = params.movDir;
            config.movInterval = params.movInterval;
            config.numMovPos = params.numMovPos;
            movStep = params.movInterval .* params.movDir./vecnorm(params.movDir,2);
            
            config.rcvPos = params.rcvPos;
            config.srcPos = ones(params.numMovPos,3).*params.srcPos +...
                (0:(params.numMovPos-1)).' * movStep;
            
        otherwise % just 1 src and 1 rcv position, nothing moving
            config.srcPos = params.srcPos;
            config.rcvPos = params.rcvPos;
            config.numMovPos = 1;
    end
    
    global globalSurfaceCount
    PCs = struct('n',cell(1,config.numMovPos),...
                'pos',cell(1,config.numMovPos),...
                'mass',cell(1,config.numMovPos));

    for ii = 1:config.numMovPos
        globalSurfaceCount = 1;
        ISM = IsmSimulator();
        ISM.fs = params.fs;
        ISM.geomFile = params.filename;
        ISM.maxIsOrder = params.maxIsOrder;
        
        switch params.moving
            case 'rcv'
                ISM.sourcePosition = config.srcPos;
                ISM.receiverPosition = config.rcvPos(ii,:);
            case 'src'
                ISM.sourcePosition = config.srcPos(ii,:);
                ISM.receiverPosition = config.rcvPos;
            otherwise
                ISM.sourcePosition = config.srcPos;
                ISM.receiverPosition = config.rcvPos;
        end
        
        PCs(ii) = ISM.simulateIs;
        
        if params.doPlot ~= 0
            if ii == 1 || ii == round(config.numMovPos/2) || ii == config.numMovPos
                ISM.plotGeometry;
                view([0 90])
            end
        end
        clear ISM;
    end
end