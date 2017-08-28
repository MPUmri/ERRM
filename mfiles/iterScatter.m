function iterScatter(x,ys,sz,colors,lineWidth)
    % Make multiple scatter plots
    % Inputs:
    % x: vector [Mx1], values for x-axis
    % ys: matrix [MxN], N sets of values for y-axis
    % colours: cell [Nx1] (optional), cell containing colour for each line
    % lineWidth: real (optional), width of the lines

    if length(x) ~= size(ys,1)
        ys = ys';
    end
    N = size(ys,2);
    
    if nargin<5
        lineWidth = 2;
    end
    if nargin<4
        for i=1:N
            colors{i} = rand(1,3);
        end
    end
    
    hold on
    for i=1:N
       scatter(x,ys(:,i),sz,colors{i},'LineWidth',lineWidth);
    end
    hold off

end

