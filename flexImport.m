function data = flexImport(fname)
% load data located in fname to variable for multiple common formats
% supported file types: *.mat, *.tiff, *.tif


[~,~,ext] = fileparts(fname);

switch ext
    
    case ".tiff"
        data = Tiff2Var(fname);
    
    case ".tif"
        data = Tiff2Var(fname);
        
    case ".mat"
        data = importdata(fname);


                %If loaded as a struct, 
        if isstruct(data)
        varName = fieldnames(data);
        data = data.(varName{1});
%         if length(varName)>1; throw(MException('myComponent:inputError', 'Imported Data has too many fields')); end
        if length(varName)>1; warning('Multiple fields detected, only loading first.');
        end
        
end

end
