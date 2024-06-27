function [] = view_flux_max(filename,color_bins)
%--------------------------------------------------------------------------
%[] = view_flux_integrated(filename,color_bins)
% Plots vectors for the time integrated reaction fluxes.  
% Up to 6 sets of fluxes are plotted in different colors.
% Inputs>  filename: file from which flux data is read.
%          cycle_limit: The limits of integration, [1 0] inclues all.
%          color_bins: array specifing the range of flux for each color. 
% Outputs: None
%--------------------------------------------------------------------------
  font_size = 16;

% Read TS file
  [zz, aa, xmf, time, temperature, density, timestep, ~, flx_end, flx] = read_ts_file(filename);
  nn=aa-zz;

% Plot nuclear chart
  [z_max,z_min,n_max,n_min,point_size] = draw_nz_background(zz,nn,font_size);

% Identify starting and ending points for all fluxes
  zt=zz(flx_end(:,1));  nt=nn(flx_end(:,1));
  zp=zz(flx_end(:,2));  np=nn(flx_end(:,2));
  
% Find max fluxes
  flx_write = zeros(size(flx,1),1);
  [~,iwrite] = max(abs(flx),[],2);
  for j = 1:size(flx,1)
      flx_write(j) = flx(j,iwrite(j));
  end
  
% Label Figure
  cond_label={['Maximum flux ',num2str(1,'%5d'),' to ',num2str(size(time,2),'%5d')],...
    ['Time =',num2str(1,'%8.3e'),' to ',num2str(size(time,2),'%8.3e'),' seconds']};
   text(n_min+1, z_max-1.5,cond_label,'FontSize',14,'HorizontalAlignment','left');

% Choose maximum flux
  flx_max=max(abs(flx_write));

% Draw flux vectors  
  draw_nz_flux(flx_max,color_bins,flx_write,zt,nt,zp,np)
  
end