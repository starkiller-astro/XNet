function [] = view_mf_max(filename,x_cut)
%--------------------------------------------------------------------------
%[] = view_mf_max(filename,x_cut)
% Scatter plot color coded by maximum mass fractions during evolution.  
% Inputs>  filename: file from which flux data is read.
%          x_cut: minimum mass fraction to plot.
% Outputs: None
%--------------------------------------------------------------------------
  font_size = 16
% Read TS file
  [zz, aa, xmf, time, temperature, density, ~, ~, ~, ~] = read_ts_file(filename);
  nn=aa-zz;

% Plot nuclear chart
  [z_max,z_min,n_max,n_min,point_size] = draw_nz_background(zz,nn,font_size);
  point_area = point_size^2;


% Find maximim mass fraction for each species
  x_max           =max(xmf,[],2);

% Label Figure
  text(n_min+1, z_max-1.5,'Maximum Mass Fraction','FontSize',14,'HorizontalAlignment','left');

% Plot points with color according to mass fraction 
  draw_nz_mf(nn,zz,x_max,x_cut,font_size,point_area)
  hold off
    
end