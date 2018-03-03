function [imgData, hdr, acqTime] = getDicomImage(dcmDir)

    dcmFiles = dir([dcmDir '/*.dcm']);

    tmpImg = dicomread(fullfile(dcmDir,dcmFiles(1).name));
    [sX, sY] = size(tmpImg);
    sZ = length(dcmFiles);

    imgData = zeros(sX,sY,sZ);

    acqTime = zeros(sZ,1);
    for i=1:length(dcmFiles)
       hdr = dicominfo(fullfile(dcmDir,dcmFiles(i).name));
       imgData(:,:,hdr.InstanceNumber) = dicomread(hdr);
       acqTime(hdr.InstanceNumber) = str2num(hdr.AcquisitionTime);
    end

end
