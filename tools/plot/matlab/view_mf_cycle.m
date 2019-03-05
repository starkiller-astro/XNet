function [] = view_mf_cycle(filename,cycle_number,x_cut)
%--------------------------------------------------------------------------
%[] = view_mf_cycle(filename,cycle_number,x_cut)
% Scatter plot color coded by mass fractions at a single cycle_number.  
% Inputs>  filename: file from which flux data is read.
%          cycle_number: Identity of cycle to be plotted
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

% Choose cycle to plot
  x_cycle           =xmf(:,cycle_number);
  time_cycle        =time(cycle_number);
  temperature_cycle =temperature(cycle_number);
  density_cycle     =density(cycle_number);

% Label Figure
  cond_label={['Time =',num2str(time_cycle,'%8.3e'),' seconds'],['T= ',...
      num2str(temperature_cycle,'%5.2f'),' GK'],[' \rho= ', num2str(density_cycle,'%8.2e'),' g/cm^3']}
  text(n_min+1, z_max-1.5,cond_label,'FontSize',16,'HorizontalAlignment','left')

% Plot points with color according to mass fraction 
  draw_nz_mf(nn,zz,x_cycle,x_cut,font_size,point_area)
  hold off
  
end