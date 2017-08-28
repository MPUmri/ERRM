function [dirLocation] = DefaultFolders()
    % This sets the default folder locations

    dirLocation.sim = './data/simData';
    dirLocation.qiba = './data/QIBA-ToftsV6';
    dirLocation.qinZip = './data/QIN Breast DCE-MRI';
    dirLocation.qin = './data/QINBreast';
    dirLocation.sarc = './data/sarcData';

    dirLocation.results = './dataResults';
end
