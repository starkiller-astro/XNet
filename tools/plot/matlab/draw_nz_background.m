function [z_max,z_min,n_max,n_min,point_size] = draw_nz_background(zz,nn,font_size)
%--------------------------------------------------------------------------
% [z_max,z_min,n_max,n_min,point_size] = draw_nz_background(zz,nn,font_size)
% Place points to mark species included in network.  Place elemental and 
% atomic mass labels.
% Inputs>  zz: proton number of included species
%          nn: neutron number of included species
%          font_size: size of font for axis labels
% Outputs> z_max: maximum proton number in network
%          z_min: minimum proton number in network 
%          n_max: maximum neutron number in network 
%          n_min: minimum neutron number in network 
%          point_size: size of points that prevents overlap
%--------------------------------------------------------------------------
% Find all included isotopes for each species.
%  zall=[zt;zp]; nall=[nt;np];
  z_max=max(zz); z_min=min(zz);
  n_max=max(nn); n_min=min(nn);
  zpoint=[]; npoint=[];

% Determine Optimal point and label sizes.
    set(gca, 'Units', 'Points');
    axis_position = get(gca, 'Position');
    axis_xlen = axis_position(3) - axis_position(1);
    axis_ylen = axis_position(4) - axis_position(2);
    zlen = z_max - z_min + 2;
    nlen = n_max - n_min + 2;
    point_size_x = axis_xlen / nlen;
    point_size_y = axis_ylen / zlen;
    point_size = min(point_size_x,point_size_y);
    label_size = fix(point_size_y);
    if (label_size > font_size) ;
      label_size = font_size;
    else if (label_size < 3); 
      label_size = 3;
      end;
    end;

% Loop over all elements, including neutrons (z=0)
  for i= 1:z_max+1
    zint(i)=i-1;
    iz=find(zz==i-1);
%   Find the lightest and heaviest isotope 
    if not(isempty(iz))
      leftn(i)=min(nn(iz));
      rightn(i)=max(nn(iz));
      alabel(i)=(rightn(i)+zint(i));
      for j=leftn(i):rightn(i)
        in=find(nn(iz)==j);
        if not(isempty(in))
          zpoint=[zpoint;i-1];
          npoint=[npoint;j];
        end
      end
    end
  end

% Plot species
  h=scatter(npoint,zpoint,10,'k','.');
  plot_id(1)=h(1);

% Overplot
  hold on

  axis([n_min-2 n_max+2 z_min-1 z_max+1]) 
  axis equal
  set(gca,'FontSize',font_size)
  nshell=[2,8,20,28,50,82,126];
  xlabel('N (Neutron Number)')
  set(gca,'xtick',nshell)
  ylabel('Z (Proton Number)')
  pshell=[2,8,20,28,50,82];
  set(gca,'ytick',pshell)
  grid on

% Load element names
  element=build_element_symbol();

% Label Elements to the left
  text(leftn-.5,zint,element(1:z_max+1,:),'FontSize',label_size,'HorizontalAlignment','right','Clipping','on')

% Label Mass numbers (for Even Z only)
  acell=num2cell(alabel);
  text(rightn(3:2:end)+.5,zint(3:2:end)-.5,acell(3:2:end),'FontSize',label_size,'HorizontalAlignment','left','Rotation',-45,'Clipping','on');

end