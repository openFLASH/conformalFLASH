%% exportCEFtoMCNPX
% Function exporting ridge filter coordonates to MCNPX text format
%
%% Syntax
% |exportCEFtoMCNPX(coordFilename, GridLayout, SpikeCentres, BaseSize, SpikeSteps, CylinderRadii)|
%
%% Description
% |exportCEFtoMCNPX(coordFilename, GridLayout, SpikeCentres, BaseSize, SpikeSteps, CylinderRadii)| Export ridge filter elemental coordinates to .i text format in a folder specified by user.
%
% * c RPP   RECTANGULAR PARALLELEPIPED
% 11  RPP  -20.00  20.00  $ xmin xmax
%          -20.00  20.00  $ ymin ymax
%          -20.00  20.00  $ zmin zmax
%
% * c Centered regular hexagon (CRH)
% 11  CRH  0.00  0.00  0.00       $ x, y, z coordinates of the bottom of the hexagonal prism
%          0.00  0.00  2.00       $ vector from the bottom to the top of the hexagonal prism
%          8.66  -5.00  0.00      $ vector to center of the 1st facet
%          0.00  -10.00  0.00     $ vector to center of the 2nd facet
%          -8.66  -5.0000  0.00   $ vector to center of the 3rd facet
%
% * c Spike centered regular hexagon/square (SCRH/SCS)
% 11-2  SCRH  0.00  0.00  2.00      $ x, y, z coordinates of the bottom of the cylinder
%             0.00  0.00  5.00      $ vector from the bottom to the top of the cylinder
%             8.66  -5.00  0.00     $ vector to center of the 1st facet
%             0.00  -10.00  0.00    $ vector to center of the 2nd facet
%             -8.66  -5.0000  0.00  $ vector to center of the 3rd facet
%% Input arguments
% |coordFilename| - _STRING_ - Location where .i MNCPX coordinate files for each beam are saved to
%
% |GridLayout| - _STRING_ - Type of spot grid currently available - SQUARE or HEXAGONAL
%
% |SpikeCentres| - _VECTOR of VECTOR_ - Vector containing the central x,y coodinates of each vertical structure in the ridge filter. The x,y coordinates are are expressed in [mm] in IEC Gantry Coordinates System
%
% |BaseSize| - _STRUCT of VECTOR_ - structure defining the characteristics of the base
%     * |BaseSize.apothem(i)| -_SCALAR VECTOR_- Apothem (mm) of the base of the i-th spike
%     * |BaseSize.baseHeight(i)| -_SCALAR VECTOR_- Maximum height (mm) of the i-th spike
%
% |SpikeSteps| - _CELL of VECTOR_ - cell containing the levels of each spike in the ridge filter structure
%
% |PolygonApothems| - _CELL of VECTOR_ - cell containing the apothems per level of each spike in the ridge filter structure
%
%% REFERENCE
%
%
%% Contributors
% Authors : Lucian Hotoiu (open.reggui@gmail.com)

function exportCEFtoMCNPX(coordFilename, GridLayout, SpikeCentres, BaseSize, SpikeSteps, PolygonApothems, spikeDirection)

    bases = {};
    spikes = {};

    if strcmp(GridLayout, 'SQUARE')
        % Compute deviation from center of square to lateral face in x, y z direction. The same for all faces and all hexagonal bases
        apothemSquare = BaseSize.apothem;

        xmin = SpikeCentres(:,1) - apothemSquare';
        xmax = SpikeCentres(:,1) + apothemSquare';
        ymin =  SpikeCentres(:,2) - apothemSquare';
        ymax =  SpikeCentres(:,2) + apothemSquare';
        zmin = zeros(length(SpikeCentres),1);
        zmax = BaseSize.baseHeight';

        % Structure coordinates per pair min/max
        x = [xmin, xmax];
        y = [ymin, ymax];
        z = [zmin, zmax];

        % Make all square bases coordinates
        bases = {x, y, z};


        %------------------------------------------------------------------
        % Pile up square towers coordinates
        for spike = 1:length(SpikeSteps)
            xs{spike} = [];
            ys{spike} = [];
            zs{spike} = [];
            sxmin = [];
            sxmax = [];
            symin = [];
            symax = [];
            szmin = [];
            szmax = [];

            z_prev = BaseSize.baseHeight(spike);

            for step = 2:length(cell2mat(SpikeSteps(spike)))
                apothemSquare = PolygonApothems{spike}(step);

                if (apothemSquare ~= 0)
                    sxmin = [sxmin; SpikeCentres(spike,1) - apothemSquare];
                    sxmax = [sxmax; SpikeCentres(spike,1) + apothemSquare];
                    symin = [symin; SpikeCentres(spike,2) - apothemSquare];
                    symax = [symax; SpikeCentres(spike,2) + apothemSquare];
                    szmin = [szmin; z_prev];
                    szmax = [szmax; SpikeSteps{spike}(step)];
                    z_prev = szmax(step-1);
                end
            end

            xs{spike} = [xs{spike},[sxmin, sxmax]];
            ys{spike} = [ys{spike},[symin, symax]];
            zs{spike} = [zs{spike},[szmin, szmax]];
        end

        % Make all square spikes coordinates
        spikes = {xs, ys, zs};

    elseif strcmp(GridLayout, 'HEXAGONAL')
        % Compute deviation from center of hexagon to lateral face. The same for all faces and all hexagonal bases
        angles = [30,60,90]; % angles in the triangle containing the vector from center to the 1st and 3rd facet
        apothemHex = BaseSize.apothem;
        heightBase = BaseSize.baseHeight;
        %side = 2.* apothemHex ./ sqrt(3);
        [dx, dy] = computeTriangle(apothemHex, angles);

        % define the origin point for the directional vectors as the central axis point of all hexagons
        x0 = SpikeCentres(:,1);
        y0 = SpikeCentres(:,2);
        z0 = zeros(size(SpikeCentres,1),1);
        O0 = [x0, y0, z0];

        % compute the directional vectors from each hexagonal central axis point to its 3 right facets
        % facet 1 - up right
        % central point vector to 1st facet is O0F1 = F1 - O0
        x1 = SpikeCentres(:,1) + dx';
        y1 = SpikeCentres(:,2) + dy';
        z1 = z0;
        F1 = [x1, y1, z1];
        O0F1 = F1 - O0;

        % facet 2 - center right
        % central point vector to 2st facet is O0F2 = F2 - O0
        x2 = SpikeCentres(:,1) + apothemHex';
        y2 = SpikeCentres(:,2);
        z2 = z0;
        F2 = [x2, y2, z2];
        O0F2 = F2 - O0;

        % facet 3 - bottom right
        % central point vector to 1st facet is O0F3 = F3 - O0
        x3 = SpikeCentres(:,1) + dx';
        y3 = SpikeCentres(:,2) - dy';
        z3 = z0;
        F3 = [x3, y3, z3];
        O0F3 = F3 - O0;

        % top face
        % central point vector to top facet is O0FT = FT - O0
        xt = SpikeCentres(:,1);
        yt = SpikeCentres(:,2);
        zt = heightBase';
        FT = [xt, yt, zt];
        O0FT = FT - O0;

        % Make all hex bases coordinates
        bases = {O0, O0FT, O0F1, O0F2, O0F3};


        %------------------------------------------------------------------
        % Pile up hexagon towers coordinates
        for spike = 1:length(SpikeSteps)
            hb{spike} = [];
            ht{spike} = [];
            vt1{spike} = [];
            vt2{spike} = [];
            vt3{spike} = [];

            for step = 2:length(cell2mat(SpikeSteps(spike)))
                % all hexagon vectors to faces
                apothemHex = PolygonApothems{spike}(step);
                %side = 2.* apothemHex ./ sqrt(3);
                [dx, dy] = computeTriangle(apothemHex, angles);

                if (apothemHex ~= 0)
                    xt0 = SpikeCentres(spike,1);
                    yt0 = SpikeCentres(spike,2);
                    zt0 = SpikeSteps{spike}(step-1);
                    Ot0 = [xt0, yt0, zt0];

                    % compute the vectors from each hexagonal centre to its 3 left facets
                    % facet 1 - up right
                    xt1 = SpikeCentres(spike,1) + dx';
                    yt1 = SpikeCentres(spike,2) + dy';
                    zt1 = zt0;
                    Ft1 = [xt1, yt1, zt1];
                    Ot0Ft1 = Ft1 - Ot0;

                    % facet 2 - center right
                    xt2 = SpikeCentres(spike,1) + apothemHex';
                    yt2 = SpikeCentres(spike,2);
                    zt2 = zt0;
                    Ft2 = [xt2, yt2, zt2];
                    Ot0Ft2 = Ft2 - Ot0;

                    % facet 3 - bottom right
                    xt3 = SpikeCentres(spike,1) + dx';
                    yt3 = SpikeCentres(spike,2) - dy';
                    zt3 = zt0;
                    Ft3 = [xt3, yt3, zt3];
                    Ot0Ft3 = Ft3 - Ot0;


                    % define vectors to center of each of the 3 left hexagon facets
                    % vector facet 1 - up right
                    vt1{spike} = [vt1{spike}; Ot0Ft1];
                    % vector facet 2 - center right
                    vt2{spike} = [vt2{spike}; Ot0Ft2];
                    % vector facet 3 - bottom right
                    vt3{spike} = [vt3{spike}; Ot0Ft3];

                    % all hexagon top vectors
                    xTt = xt0;
                    yTt = yt0;
                    zTt = SpikeSteps{spike}(step);
                    FTt = [xTt, yTt, zTt];
                    O0FTt = FTt - Ot0;

                    ht{spike} = [ht{spike}; O0FTt];
                end
            end
            % all hexagon base coordinates
            if (~isempty(ht{spike}))
                hb{spike} = ht{spike};
                hb{spike}(end,:) = [];
                % add hexagonal-base top coordinate as base coordinate of the 1st peak in spike
                hb{spike} = [[SpikeCentres(spike,:), BaseSize.baseHeight(spike)] ; hb{spike}];
            end
        end

        % Make all hex spikes coordinates
        spikes = {hb, ht, vt1, vt2, vt3};
    end

    % write the text file in MNCPX format
    writeMNCPXFormat(coordFilename, GridLayout, bases, spikes);
end

function [x,y] = computeTriangle(apothem, angles)
    % The small triangle delimited in between the hexagon corner and the
    % vertical line down from half length its side
    A = angles(1);
    B = angles(2);
    C = angles(3);
    x = (sin(A*pi/180).*apothem)/sin(C*pi/180);
    y = (sin(B*pi/180).*apothem)/sin(C*pi/180);
end


function writeMNCPXFormat(coordFilename, GridLayout, bases, spikes)

    fid = fopen(coordFilename,'wt');

    % write bases
    if strcmp(GridLayout, 'SQUARE')
        for iBase = 1:size(bases{1},1)
            % write xmin xmax coordinates
            formatSpec = 'c \nc RECTANGULAR PARALLELEPIPED \n%s  RPP  %3.9f  %3.9f\n';
            fprintf(fid,formatSpec,num2str(iBase),bases{1}(iBase,:));

            % write ymin ymax coordinates
            formatSpec = '           %3.9f  %3.9f\n';
            fprintf(fid,formatSpec,bases{2}(iBase,:));

            % write zmin zmax coordinates
            formatSpec = '           %3.9f  %3.9f\n';
            fprintf(fid,formatSpec,bases{3}(iBase,:));

        end
        %--------------------------------------------------------------

        % Write square spikes
        for iSpike = 1:length(spikes{1})
            for jLevel = 1:size(spikes{1,1}{iSpike},1)
                % write xmin xmax coordinates
                formatSpec = 'c \nc Spike RECTANGULAR PARALLELEPIPED \n%s-%s SCS  %3.9f  %3.9f\n';
                fprintf(fid,formatSpec,num2str(iSpike), num2str(jLevel), spikes{1,1}{iSpike}(jLevel,:));

                % write ymin ymax coordinates
                formatSpec = '           %3.9f  %3.9f\n';
                fprintf(fid,formatSpec,spikes{1,2}{iSpike}(jLevel,:));

                % write zmin zmax coordinates
                formatSpec = '           %3.9f  %3.9f\n';
                fprintf(fid,formatSpec,spikes{1,3}{iSpike}(jLevel,:));
            end
        end
    elseif strcmp(GridLayout, 'HEXAGONAL')
        % Write hexagonal bases
        for iBase = 1:size(bases{1},1)
            % write bottom face coordinates
            formatSpec = 'c \nc Centered regular hexagon \n%s  CRH  %3.9f  %3.9f  %3.9f\n';
            fprintf(fid,formatSpec,num2str(iBase),bases{1}(iBase,:));

            % write top face coordinates
            formatSpec = '           %3.9f  %3.9f  %3.9f\n';
            fprintf(fid,formatSpec,bases{2}(iBase,:));

            % write vector to facet 1 - upper left
            formatSpec = '           %3.9f  %3.9f  %3.9f\n';
            fprintf(fid,formatSpec,bases{3}(iBase,:));

            % % write vector to facet 2 - left
            formatSpec = '           %3.9f  %3.9f  %3.9f\n';
            fprintf(fid,formatSpec,bases{4}(iBase,:));

            % % write vector to facet 3 - lower left
            formatSpec = '           %3.9f  %3.9f  %3.9f\n';
            fprintf(fid,formatSpec,bases{5}(iBase,:));
        end


        %------------------------------------------------------------------

        % Write hexagonal spikes
        for iSpike = 1:length(spikes{1})
            for jLevel = 1:size(spikes{1,1}{iSpike},1)
                % write bottom face coordinates
                formatSpec = 'c \nc Spike centered regular hexagon \n%s-%s SCRH  %3.9f  %3.9f  %3.9f\n';
                fprintf(fid,formatSpec,num2str(iSpike), num2str(jLevel),spikes{1,1}{iSpike}(jLevel,:));

                % write top face coordinates
                formatSpec = '           %3.9f  %3.9f  %3.9f\n';
                fprintf(fid,formatSpec,spikes{1,2}{iSpike}(jLevel,:));

                % write vector to facet 1 - upper left
                formatSpec = '           %3.9f  %3.9f  %3.9f\n';
                fprintf(fid,formatSpec,spikes{1,3}{iSpike}(jLevel,:));

                % % write vector to facet 2 - left
                formatSpec = '           %3.9f  %3.9f  %3.9f\n';
                fprintf(fid,formatSpec,spikes{1,4}{iSpike}(jLevel,:));

                % % write vector to facet 3 - lower left
                formatSpec = '           %3.9f  %3.9f  %3.9f\n';
                fprintf(fid,formatSpec,spikes{1,5}{iSpike}(jLevel,:));

                % write radii
                %formatSpec = '           %3.9f\n';
                %fprintf(fid,formatSpec,spikes{1,3}{iSpike}(jLevel));
            end
        end
    end

    fclose(fid);
end
