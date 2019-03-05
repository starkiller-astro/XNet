function [] = plot_z_mf(filename,varargin)
%--------------------------------------------------------------------------
%[] = plot_z_mf(filename,MassLimit,CycleNumber)
% Scatter plot color coded by mass fractions at a single cycle_number.  
% Inputs>  filename: file from which data is read.
%          MassLimit: lower limit of mass fractions to plot.
%          CycleNumber: Identity of cycle to be plotted, default is final
% Outputs: None
%--------------------------------------------------------------------------

% Create an instance of the inputParser class.
  p = inputParser;

% Define required inputs
  p.addRequired('filename',@ischar);

% Define optional inputs
  p.addOptional(  'MassLimit', 1e-25, @(x)validateattributes(x, {'numeric'}, ...
    {'scalar', 'real', 'positive', '>=', 1e-30, '<=', 1e-1}));
  p.addOptional('CycleNumber', 0, @(x)validateattributes(x, {'numeric'}, ...
    {'scalar', 'integer'}));

% Parse and validate all input arguments.
  p.parse(filename, varargin{:});

% Assign input arguments
  filename      = p.Results.filename;
  mass_limit    = p.Results.MassLimit;
  cycle_number  = p.Results.CycleNumber;

% Read TS file
  [zz, aa, xmf, time, temperature, density, ~, ~, ~, ~] = read_ts_file(filename);
  nn=aa-zz;
  if (cycle_number == 0) 
    cycle_number=size(time,2)
  end
  
% Choose cycle to plot
  x_cycle           =xmf(:,cycle_number);
  time_cycle        =time(cycle_number);
  temperature_cycle =temperature(cycle_number);
  density_cycle     =density(cycle_number);

% Sum mass fractions over Z,
  zint=cast(zz,'int32');
  zint=zint+1; % Z=0 is possible, so shift Z array to right by one
  zmin=min(zint);
  zmax=max(zint); 
  xz = zeros(1,zmax);
  ny = size(x_cycle);
  for n =1:ny
    xz(zint(n))=xz(zint(n))+x_cycle(n);
  end

% Plot summed mass fractions, shifted to left be one to accomodate Z=0. 
  zmin=zmin-1; 
  zmax=zmax-1;
  semilogy(zmin:zmax,xz);
  ylim([mass_limit 1.1])
  xlim([zmin zmax])
  xlabel('Z (Proton Number)')
  ylabel('Mass Fraction')
  hold on

% Label Figure
  xmax=max(xz)
  cond_label={['Time =',num2str(time_cycle,'%8.3e'),' seconds'],['T= ',...
      num2str(temperature_cycle,'%5.2f'),' GK'],[' \rho= ', num2str(density_cycle,'%8.2e'),' g/cm^3']}
  text(5,0.1,cond_label,'FontSize',14,'HorizontalAlignment','left')
  
  hold off
  
end