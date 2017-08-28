function shadedErrorPlot(x,y,lowErr,hiErr,c,doScatter)

    if nargin<6
        doScatter = false;
    end

    if nargin<5
        c = [.9 .9 .9];
    end

    %%
    if doScatter
       scatter(x,y,'filled','MarkerEdgeColor',c,'MarkerFaceColor',c); 
    end
    %%
    fill([x;flipud(x)],[lowErr;flipud(hiErr)],c,'linestyle','none','FaceAlpha',.3);
    line(x,y,'Color',c,'LineWidth',2)
    
end

