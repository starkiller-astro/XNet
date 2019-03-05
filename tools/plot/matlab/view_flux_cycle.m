function [] = view_flux_cycle(filename,cycle_number,color_bins)
%--------------------------------------------------------------------------
%[] = view_flux_cycle(filename,cycle_number,color_bins)
% Plots vectors for the reaction fluxes at a single cycle_number.  
% Up to 6 sets of fluxes are plotted in different colors.
% Inputs>  filename: file from which flux data is read.
%          cycle_number: Identity of cycle to be plotted
%          color_bins: array specifing the range of flux for each color. 
% Outputs: None
%--------------------------------------------------------------------------
  font_size = 16

% Read TS file
  [zz, aa, xmf, time, temperature, density, ~, ~, flx_end, flx] = read_ts_file(filename);
  nn=aa-zz;

% Plot nuclear chart
[z_max,z_min,n_max,n_min,point_size] = draw_nz_background(zz,nn,font_size);

% Identify starting and ending points for all fluxes
  zt=zz(flx_end(:,1));  nt=nn(flx_end(:,1));
  zp=zz(flx_end(:,2));  np=nn(flx_end(:,2));
  
% Choose cycle to plot
  flx_write         =flx(:,cycle_number)
  time_write        =time(cycle_number);
  temperature_write =temperature(cycle_number);
  density_write     =density(cycle_number);
  
% Label Figure
  cond_label={['Time =',num2str(time_write,'%8.3e'),' seconds'],['T= ',...
      num2str(temperature_write,'%5.2f'),' GK'],[' \rho= ', num2str(density_write,'%8.2e'),' g/cm^3']}
  text(n_min+1, z_max-1.5,cond_label,'FontSize',14,'HorizontalAlignment','left')

% Choose maximum flux
  flx_max=max(abs(flx_write))

% Draw flux vectors  
  draw_nz_flux(flx_max,color_bins,flx_write,zt,nt,zp,np)
  hold off
  
end