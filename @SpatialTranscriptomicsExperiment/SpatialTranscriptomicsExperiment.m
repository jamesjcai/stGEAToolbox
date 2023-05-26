classdef SpatialTranscriptomicsExperiment
    properties
        sce SingleCellExperiment   %single-spot experiment is an SCE
        img                        % Image data as an array img=imread(image_paths);
        xy                         % initial position of spot spatial map
        tissue_positions_list
        scalefactors_json          % https://support.10xgenomics.com/spatial-gene-expression/software/pipelines/latest/output/images
    end
    properties (Dependent)
      NumSpots
      NumGenes
    end

   methods
      function obj = SpatialTranscriptomicsExperiment(sce,xy,img)
          if nargin<1 || isempty(sce), sce=SingleCellExperiment(); end          
          if nargin<2 || isempty(xy), xy=[rand(sce.NumCells,1) rand(sce.NumCells,1)]; end
          if nargin<3 || isempty(img), img=[]; end
          %if nargin<4 || isempty(scalef), scalef=1.0; end
          %obj.scalef=scalef;
          if ~any(strcmp(sce.list_cell_attributes,'ste.xy'))
            sce.list_cell_attributes{end+1}='ste.xy';
            sce.list_cell_attributes{end+1}=xy;
          end
          obj.sce=sce;
          obj.xy=xy;
          obj.img=img;          
      end

%    function out = get.s(obj)
%       out = obj.scalef*((obj.sce.s-min(obj.sce.s))./(max(obj.sce.s)-min(obj.sce.s)));
%       out = out - obj.xy;
%    end

    obj = writepositions(obj,filename)

    function obj = rescalexy(obj,scalef)
        obj.xy=obj.xy*scalef;
        %if nargin<2, scalef=1.0; end
        %obj.xy=scalef*((obj.xy-min(obj.xy))./(max(obj.xy)-min(obj.xy)));        
    end

    function obj = rotateexy(obj,theta)
        if nargin<2, theta=0.5*pi; end
        a=obj.xy-mean(obj.xy);
        obj.xy=([cos(theta) -sin(theta); sin(theta) cos(theta)]*a')'+mean(obj.xy);
    end

    function obj = zoomxy(obj,f)
        if nargin<2, f=0.05; end
        a=obj.xy-mean(obj.xy);
        a=(1+f)*a;
        obj.xy=a+mean(obj.xy);      
    end    

    function r = title(obj)
       r=sprintf('%d x %d\n[genes x spots]',...
           size(obj.sce.X,1),size(obj.sce.X,2));
    end

   function m = get.NumSpots(obj)
      m = size(obj.sce.X,2); 
   end
   
   function obj = set.NumSpots(obj,~)
      fprintf('%s%d\n','NumSpots is: ',obj.NumSpots)
      error('You cannot set NumSpots property'); 
   end
   
   function m = get.NumGenes(obj)
      m = size(obj.sce.X,1); 
   end
   
   function obj = set.NumGenes(obj,~)
      fprintf('%s%d\n','NumGenes is: ',obj.sce.NumGenes)
      error('You cannot set NumGenes property'); 
   end
 
    function r=numspots(obj)
        r=size(obj.sce.X,2);
    end
    
    function r=numgenes(obj)
        r=size(obj.sce.X,1);
    end

   end
end