function iterPlot(x,ys,colors,lineStyle,lineWidth)
    % Make multiple plots
    % Inputs:
    % x: vector [Mx1], values for x-axis
    % ys: matrix [MxN], N sets of values for y-axis
    % colours: cell [Nx1] (optional), cell containing colour for each line
    % lineStyle: cell [Nx1] (optional), cell containing line style for each line
    % lineWidth: real (optional), width of the lines

    if length(x) ~= size(ys,1)
        ys = ys';
    end
    N = size(ys,2);
    
    if nargin<5
        lineWidth = 2;
    end
    if nargin<4
        lineStyleList = ['-','--',':','-.'];
        for i=1:N
            lineStyle{i} = '-';
            % If random lineStyles desired, use:
            % lineStyle{i} = randsample(lineStyleList,1);
        end
    end
    if nargin<3
        for i=1:N
            colors{i} = rand(1,3);
        end
    end
    
    hold on
    for i=1:N
       plot(x,ys(:,i),lineStyle{i},'Color', colors{i},'LineWidth',lineWidth);
    end
    hold off

end

