function [foundDirs] = findDir(startDir, targetString, foundDirs)
    
    if nargin<3
        foundDirs = {};
    end

    dirList = dir(startDir);
    
    targetList = find(contains({dirList.name},targetString));
    if length(targetList)
        foundDirs{end+1}=startDir;
        return
    end
    
    for i=3:length(dirList)
        foundDirs = findDir(fullfile(startDir,dirList(i).name), targetString, foundDirs);
    end
end