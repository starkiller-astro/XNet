function [axis_h,patch_h] = draw_nz_network_patch(data_directory)
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
  [minn,~] = min(nn(zz==minz));
  [maxn,~] = max(nn(zz==minz));
  xvertices = maxn+0.5:-1:minn-0.5;
  yvertices = repmat(minz-0.5,1,maxn-minn+2);
  for iz=minz:maxz
      [minn,~] = min(nn(zz==iz));
      xvertices = [xvertices,minn-0.5,minn-0.5];
      yvertices = [yvertices,iz-0.5,iz+0.5];
  end
  [minn,~] = min(nn(zz==maxz));
  [maxn,~] = max(nn(zz==maxz));
  xvertices = [xvertices,minn-0.5:maxn+0.5];
  yvertices = [yvertices,repmat(maxz+0.5,1,maxn-minn+2)];
  for iz=maxz:-1:minz
      [maxn,~] = max(nn(zz==iz));
      xvertices = [xvertices,maxn+0.5,maxn+0.5];
      yvertices = [yvertices,iz+0.5,iz-0.5];
  end
  patch_h = patch('XData',xvertices,'YData',yvertices,'FaceColor','none');
  axis_h = gca;
%   scatter(nn,zz,psize,'Marker','s','MarkerEdgeColor',pcolor,'MarkerFaceColor','none');

end