function [ croppedImage ] = AutoCrop( image, mask )
    % Crops 'image' using 'mask'
    % Any empty rows and column are removed

    if nargin<2
        mask = image;
    end
    
    if length(size(mask))==3
        mask = squeeze(sum(mask,3));
    end
    
    if length(size(image))==2
        image(~any(mask,2),:) = [];
        image(:,~any(mask,1)) = [];
    elseif length(size(image))==3
        image(~any(mask,2),:,:) = [];
        image(:,~any(mask,1),:) = [];
    else
        disp('Not supported');
        image=0;
    end
    
    croppedImage = image;

end