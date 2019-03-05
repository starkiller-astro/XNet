function [] = draw_nz_flux(flx_max,color_bins,flx,zt,nt,zp,np)
%--------------------------------------------------------------------------
% [] = draw_nz_flux(flx_max,color_bins,flx,zt,nt,zp,np)
% Draw vectors for the reaction flux on the NZ plane from reaction targut
% to product colored by relative flux intensity.
% Inputs>  flx_max: normalization for displayed fluxes.
%          color_bins: array specifing the range of flux for each color.
%          flx: reaction flux values
%          zt: proton number of flux target species (vector origin).
%          nt: neutron number of flux target species (vector origin).
%          zp: proton number of flux product species (vector terminus).
%          np: neutron number of flux product species (vector terminus).
% Outputs: none
%--------------------------------------------------------------------------
% set number of colors/levels
  colors=['k','r','g','b','c','m','y'];

% set number of arrays
  ncolors=size(color_bins);
  bin_upper=[1,color_bins];
  bin_lower=color_bins;

% Build Legend array
  legend_array=cell(1,ncolors(2)+1);

% Plot species
  h=scatter([nt;np],[zt;zp],10,'k','s');
  plot_id(1)=h(1);
  legend_array(1)={'Reaction Fluxes'};

% Calculate vector origin and lengths
  aflx=abs(flx);
  zo=zt;
  no=nt;
  zv=zp-zt;
  nv=np-nt;

% Reverse arrows for negative flux
  ineg=find(flx<0.0);
  zo(ineg)=zp(ineg); 
  no(ineg)=np(ineg);
  zv(ineg)=-zv(ineg);
  nv(ineg)=-nv(ineg);

% Loop over colors
for i=1:ncolors(2)
% Plot Vectors
  flx_top=flx_max*bin_upper(i);
  flx_bot=flx_max*bin_lower(i);
  ii=find(aflx<=flx_top & aflx>flx_bot);
  scale =0.0; % scale = 0 prevents autoscaling
  if i == 1
      lwidth= 2;
  elseif i == 2
      lwidth= 1.5;
  else
      lwidth= 1;
  end
  h=quiver(no(ii),zo(ii),nv(ii),zv(ii),scale,colors(i),'LineWidth',lwidth);
  plot_id(i+1)=h(1);
  label=['  > ',num2str(bin_lower(i),'%6.1e'),' of max'];
  legend_array(i+1)={label};
end

legend(plot_id,char(legend_array),'Location',[.50 .20 .1 .2]);
