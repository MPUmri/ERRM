function makeGridMaps(mapCell, chosenSlices, nameCell, climCell, cMap)

    if nargin < 5
        cMap = 'jet';
    end
    
    if nargin < 4
        % Find the global max/min and use it for all images
        % (I've never tested if this actually works - so it probably doesn't)
        curMax = 0;
        curMin = Inf;
        for i=1:length(mapCell)
           curData = mapCell{i};
           if max(curData(:)) > curMax
               curMax = max(curData(:));
           end
           if min(curData(:)) < curMin
               curMin = min(curData(:));
           end
        end
        for i=1:length(mapCell)
            climCell{i} = [curMax curMin];
        end
    end
    
    numRow = length(mapCell);
    numCol = length(chosenSlices);
    
    kounter = 1;
    figure
    for i=1:numRow
        curData = mapCell{i};
        for j=1:numCol
            subplot(4,3,kounter)
            imagesc(curData(:,:,chosenSlices(j)),climCell{i});
            colormap(cMap);
            axis square
            
            if j==1
                ylabel(nameCell{i});
            end
            
            kounter = kounter + 1;
        end
    end

end
