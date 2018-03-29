function [ outData ] = medfiltTemporal( imgData )

    outData = imgData;
    
    for i=1:size(imgData,4)
       outData(:,:,:,i) = medfilt3(imgData(:,:,:,i)); 
    end
    
end

