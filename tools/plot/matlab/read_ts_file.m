function [zz, aa, xmf, time, temperature, density, timestep, edot, flx_end, flx] = read_ts_file( ts_filename )
%--------------------------------------------------------------------------
%[zz, aa, xmf, time, temperature, density, timestep, edot, flx_end, flx] = read_ts_file( ts_filename ) 
% Reads XNet binary output file.
% Inputs>  ts_filename: name of binary file 
% Outputs< zz: proton number for each isotope in the network.
%          aa: total nucleon number for each isotope in the network
%          xmf: time series of mass fractions 
%          time: discrete temporal evolution
%          temperature:  time series of temperature
%          density: time series of density
%          timestep:  time series of discrete timesteps 
%          edot: time series of energy generation 
%          flux_end, the starting and ending points for each reaction flux.
%          flux: time seris of the reaction fluxes in the network.
%--------------------------------------------------------------------------

  file_id = fopen(ts_filename,'rb');

% Read Run Descriptions
  record_length1=fread(file_id,1,'int32');
  desc1    =setstr(fread(file_id,80,'uchar'));
  desc2    =setstr(fread(file_id,80,'uchar'));
  desc3    =setstr(fread(file_id,80,'uchar'));
  data_desc=setstr(fread(file_id,80,'uchar'));
  record_length2=fread(file_id,1,'int32');

% Read Run Settings
  record_length1=fread(file_id,1,'int32');
  kstmx    =fread(file_id,1,'int32');
  kitmx    =fread(file_id,1,'int32');
  iweak    =fread(file_id,1,'int32');
  iscrn    =fread(file_id,1,'int32');
  iconvc   =fread(file_id,1,'int32');
  changemx =fread(file_id,1,'float64');
  tolm     =fread(file_id,1,'float64');
  tolc     =fread(file_id,1,'float64');
  yacc     =fread(file_id,1,'float64');
  ymin     =fread(file_id,1,'float64');
  tdel_mm  =fread(file_id,1,'float64');
  record_length2=fread(file_id,1,'int32');

% Read Abundance Info
  record_length1=fread(file_id,1,'int32');
  abund_file =setstr(fread(file_id,80,'uchar'));
  abund_desc =setstr(fread(file_id,80,'uchar'));
  record_length2=fread(file_id,1,'int32');

% Read Thermodynamic Info
  record_length1=fread(file_id,1,'int32');
  thermo_file =setstr(fread(file_id,80,'uchar'));
  thermo_desc =setstr(fread(file_id,80,'uchar'));
  record_length2=fread(file_id,1,'int32');

% Read Nuclear Info
  record_length1=fread(file_id,1,'int32');
  ny            =fread(file_id,1,'int32');
  zz            =fread(file_id,ny,'float64');
  aa            =fread(file_id,ny,'float64');
  record_length2=fread(file_id,1,'int32');

% Read Flux Info
  record_length1=fread(file_id,1,'int32');
  nflx            =fread(file_id,1,'int32');
  if nflx>0
    flx_end(:,1)    =fread(file_id,nflx,'int32');
    flx_end(:,2)    =fread(file_id,nflx,'int32');
  end
  record_length2=fread(file_id,1,'int32');

% Read data from each timestep
  for k= 1:kstmx
    record_length1 =fread(file_id,1,'int32');
  % If end of file, exit 
    if isempty(record_length1)
      break
    end
  % Otherwise read data
    kstep          =fread(file_id,1,'int32');
    time(k)        =fread(file_id,1,'float64');
    temperature(k) =fread(file_id,1,'float64');
    density(k)     =fread(file_id,1,'float64');
    timestep(k)    =fread(file_id,1,'float64');
    edot(k)        =fread(file_id,1,'float64');
    xmf(:,k)       =fread(file_id,ny,'float64');
    if nflx>0 
      flx(:,k)       =fread(file_id,nflx,'float64');
    end
    record_length2 =fread(file_id,1,'int32');
  end

% Convert abundance to mass fraction
  for n=1:ny
    xmf(n,:) = aa(n).*xmf(n,:);
  end

% Calculate # of timesteps
  nstep = size(time,2);
  
  fclose(file_id);
end

