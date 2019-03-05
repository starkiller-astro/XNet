function [] = draw_nz_mf(nn,zz,x_plot,x_cut,font_size,point_area)
%--------------------------------------------------------------------------
%[] = flux_vect_cycle(filename,cycle_number,color_bins) 
% Scatter plot color coded by mass fractions at a single cycle_number.  
% Inputs>  filename: file from which flux data is read.
%          cycle_number: Identity of cycle to be plotted
%          xmin: minimum mass fraction to plot. 
% Outputs: None
%--------------------------------------------------------------------------

% Logarithmic mass fraction
  lx=log10(x_plot);
  
% Limit dynamic range
  lxmin=log10(x_cut);
  icut=find(lx>lxmin);

% Plot points with color according to mass fraction 
  scatter(nn(icut),zz(icut),point_area,lx(icut),'s','filled');
  cb=colorbar('horiz','Position',[.40 .20 .35 .05],'FontSize',font_size);
  caxis([lxmin 0])
  text(.42,.25,'Log (Mass Fraction)','Units','normalized','FontSize',font_size);

end