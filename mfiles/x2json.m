function [ jsonX ] = x2json( x, jsonName )

    jsonX = jsonencode(x);
    
    if nargin>1
        jsonX = ['"' jsonName '": ' jsonX];
    end

end

