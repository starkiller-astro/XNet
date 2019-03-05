function [] = movie_mf(filename,x_cut,output_mode)
%--------------------------------------------------------------------------
%[] = flux_vect_cycle(filename,cycle_number,color_bins) 
% Animation of scatter plot color coded by mass fractions.  
% Inputs>  filename: file from which flux data is read.
%          xmin: minimum mass fraction to plot.
%          output_mode: =0, output MPEG4, =1 output series of PNG
% Outputs: None
%--------------------------------------------------------------------------
  font_size = 16;

% Read TS file
  [zz, aa, xmf, time, temperature, density, ~, ~, ~, ~] = read_ts_file(filename);
  nn=aa-zz;
  cycle_max = size(time,2)
  
% Choose output format
  if (output_mode == 0)
  
% Open movie file
    movie_filename = ['./',filename,'_movie.mp4']
    WriterObj = VideoWriter(movie_filename,'MPEG-4');
    open(WriterObj)

  else  
% Create frame image destination
    directory = ['./',filename,'_frames']
    mkdir(directory);

  end  

% Set figure size
  figure
  set(gcf,'Units','pixels','Position',[200 1000 1024 768])
  
% Loop over cycles
  for cycle_number = 1:cycle_max

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
        num2str(temperature_cycle,'%5.2f'),' GK'],[' \rho= ', num2str(density_cycle,'%8.2e'),' g/cm^3']};
    text(n_min+1, z_max-1.5,cond_label,'FontSize',14,'HorizontalAlignment','left');

% Plot points with color according to mass fraction 
    draw_nz_mf(nn,zz,x_cycle,x_cut,font_size,point_area)  

% Output frame
    if (output_mode == 0)
     frame(cycle_number) = getframe(gcf);
     writeVideo(WriterObj,frame(cycle_number))
    else
      frame_filename = [directory,'/frame',num2str(cycle_number,'%05i'),'.png'];
      print('-dpng',frame_filename);
    end
    hold off;
  
  end
  
  if (output_mode == 0)
   close(WriterObj)
  end
  %movie(frame,1);
  
end