function [ny, nnz, ridx, cidx, pb, nrext, ns1, ns2, ns3, ns4 ] = read_sparse_ind( data_directory )
%--------------------------------------------------------------------------
%[zz, aa, xmf, time, temperature, density, timestep, edot, flx_end, flx] = read_ts_file( sparse_filename ) 
% Reads XNet binary sparsity pattern file
% Inputs>  sparse_filename: name of binary file 
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

  [~,aa,~,~,~,~] = read_nz_network(data_directory);
  ny = length(aa);
 
  sparse_filename = strcat(data_directory,'/sparse_ind');
  file_id = fopen(sparse_filename,'rb');

% Read Run Descriptions
  record_length1=fread(file_id,1,'int32');
  nnz    =fread(file_id,1,'int32');
  record_length2=fread(file_id,1,'int32');
  
  record_length1=fread(file_id,1,'int32');
  ridx   =fread(file_id,nnz,'int32');
  cidx   =fread(file_id,nnz,'int32');
  pb     =fread(file_id,ny+1,'int32');
  record_length2=fread(file_id,1,'int32');
  
  record_length1=fread(file_id,1,'int32');
  nrext  =fread(file_id,4,'int32');
  record_length2=fread(file_id,1,'int32');
  
  ns1 = zeros(nrext(1),1);
  ns2 = zeros(nrext(2),2);
  ns3 = zeros(nrext(3),3);
  ns4 = zeros(nrext(4),4);
  
  record_length1=fread(file_id,1,'int32');
  ns1(:,1) =fread(file_id,nrext(1),'int32');
  ns2(:,1) =fread(file_id,nrext(2),'int32');
  ns2(:,2) =fread(file_id,nrext(2),'int32');
  record_length2=fread(file_id,1,'int32');
  
  for i = 1:3
      record_length1=fread(file_id,1,'int32');
      ns3(:,i) =fread(file_id,nrext(3),'int32');
      record_length2=fread(file_id,1,'int32');
  end
  
  for i = 1:4
      record_length1=fread(file_id,1,'int32');
      ns4(:,i) =fread(file_id,nrext(4),'int32');
      record_length2=fread(file_id,1,'int32');
  end
  
  fclose(file_id);
  
end

