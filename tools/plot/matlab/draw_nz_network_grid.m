function [axis_h,group_h,line_h] = draw_nz_network_grid(data_directory)
%--------------------------------------------------------------------------
% [axis_id] = draw_nz_network(data_directory,pcolor,psize)
% Scatter plot of the isotopes in a specified netwinv file.
% Inputs>  data_directory: directory hosting netwinv file to be plotted.
%          pcolor: point color
%          psize: point size
% Outputs: axis_id: handle of current axis
%--------------------------------------------------------------------------
% Open file
  filename=strcat(data_directory,'/netwinv');
  fileID=fopen(filename);
  
% Read network size
  dataread=textscan(fileID,'%d',1);
  ny=cell2mat(dataread);

% Skip thermodata
  dataread=textscan(fileID,'%d',1);
  
% Read nuclear names
  dataread = textscan(fileID,'%s',ny);
  nname=dataread{1};

% Read nuclear data, skipping partition function  
  file_form='%*s %f %f %f %f %f';
  for i=1:ny        
    dataread=textscan(fileID,file_form,1);
    [aa(i),zz(i),nn(i),sp(i),be(i)]=dataread{1:5};
    
    dataread=textscan(fileID,'%f %f %f %f %f %f %f %f',3);

  end
  [minz,~] = min(zz);
  [maxz,~] = max(zz);
  line_h = [];
  for iz=minz:maxz
      [minn,~] = min(nn(zz==iz));
      [maxn,~] = max(nn(zz==iz));
      if length(minn) > 0 & length(maxn) > 0
          line_h = [line_h; plot( [minn-0.5,maxn+0.5], [iz-0.5, iz-0.5], 'k' )];
          line_h = [line_h; plot( [minn-0.5,maxn+0.5], [iz+0.5, iz+0.5], 'k' )];
      end
  end
  [minn,~] = min(nn);
  [maxn,~] = max(nn);
  for in=minn:maxn
      [minz,~] = min(zz(nn==in));
      [maxz,~] = max(zz(nn==in));
      if length(minz) > 0 & length(maxz) > 0
          line_h = [line_h; plot( [in-0.5,in-0.5], [minz-0.5, maxz+0.5], 'k' )];
          line_h = [line_h; plot( [in+0.5,in+0.5], [minz-0.5, maxz+0.5], 'k' )];
      end
  end
  group_h = hggroup( 'Parent', gca );
  set(get(get(group_h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
  set(line_h,'Parent',group_h);
  axis_h = gca;

end