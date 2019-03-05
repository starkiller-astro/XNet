function [] = plot_a_mf(filename,varargin)
%--------------------------------------------------------------------------
%[] = plot_a_mf(filename,MassLimit,CycleNumber)
% Scatter plot color coded by mass fractions at a single cycle_number.  
% Inputs>  filename: file from which data is read.
%          MassLimit: lower limit of mass fractions to plot.
%          cycle_number: Identity of cycle to be plotted, default is final
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
    cycle_number=size(time,2);
  end
  
% Choose cycle to plot
  x_cycle           =xmf(:,cycle_number);
  time_cycle        =time(cycle_number);
  temperature_cycle =temperature(cycle_number);
  density_cycle     =density(cycle_number);

% Sum mass fractions over A
  aint=cast(aa,'int32');
  amax=max(aint)
  xa = zeros(1,amax);
  ny = size(x_cycle);
  for n =1:ny
    xa(aint(n))=xa(aint(n))+x_cycle(n);
  end

% Plot summed mass fractions
  semilogy(xa)
  ylim([mass_limit 1.1])
  xlim([1 amax])
  xlabel('A (Mass Number)')
  ylabel('Mass Fraction')
  hold on

% Label Figure
  xmax=max(xa)
  cond_label={['Time =',num2str(time_cycle,'%8.3e'),' seconds'],['T= ',...
      num2str(temperature_cycle,'%5.2f'),' GK'],[' \rho= ', num2str(density_cycle,'%8.2e'),' g/cm^3']}
  text(5,0.1,cond_label,'FontSize',14,'HorizontalAlignment','left')
  
  hold off
  
end