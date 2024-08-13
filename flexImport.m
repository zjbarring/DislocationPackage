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
        
end

end