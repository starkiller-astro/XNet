function [ nuc_name, xmf, time, temperature, density, timestep, edot] = read_ev_file ( ev_filename )
%--------------------------------------------------------------------------
%[ nuc_name, xmf, time, temperature, density, timestep, edot] = read_ev_file ( ev_filename ) 
% Reads XNet ASCII output file.
% Inputs>  ev_filename: name of ASCII file 
% Outputs< nuc_name: symbol for each isotope in the file.
%          xmf: time series of mass fractions 
%          time: discrete temporal evolution
%          temperature:  time series of temperature
%          density: time series of density
%          timestep:  time series of discrete timesteps 
%          edot: time series of energy generation 
%--------------------------------------------------------------------------

% Set number of species
  ny=14

% Report the absence of flux data
  nflx=0

% Read Header
  label=textread(ev_filename,'%s',ny+7);

% Read data from file as block
  raw_data=textread(ev_filename,'%f','headerlines', 1);

% Reshape data into columns
  data_size=size(raw_data);
  read_width=ny+8;
  read_length=data_size(1)/read_width;
  read_data=reshape(raw_data,read_width,read_length);

% Extract variables from columns
  k=read_data(1,:);
  time=read_data(2,:);
  temperature=read_data(3,:);
  density=read_data(4,:);
  edot=read_data(5,:);
  timestep=read_data(6,:);
  xmf=read_data(7:6+ny,:);

% Reformat isotope symbols
  zz=regexprep(label(7:6+ny),'\d+','');
  aa=regexprep(label(7:6+ny),'\D+','');
  zz=regexprep(zz,'(^.)','${upper($1)}');
  nuc_name=strcat('^{',aa,'}',zz);

% Calculate # of timesteps
nstep = size(time,2)

end

