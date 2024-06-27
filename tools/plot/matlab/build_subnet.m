function [zn_new,nname_new] = build_subnet(filename,min_flx,min_xmf,max_z)

[zz, aa, xmf, time, ~, ~, timestep, ~, flx_end, flx] = read_ts_file(filename);
cycle_limit = [1, size(time,2)];

nn = aa - zz;
zn = [zz, nn];

% Maximum Z
zn_zmax_new = zn(zz <= max_z,:);

% Maximum of each species during integration
xmf_max = max(xmf,[],2);
zn_xmf_new = zn(xmf_max > min_xmf,:);

% Max flux
flx_max = zeros(size(flx,1),1);
[~,iwrite] = max(abs(flx),[],2);
for j = 1:size(flx,1)
    flx_max(j) = flx(j,iwrite(j));
end
iflx_new = find(abs(flx_max) > min_flx*max(abs(flx_max)));
znt = [zz(flx_end(:,1)), nn(flx_end(:,1))];
znp = [zz(flx_end(:,2)), nn(flx_end(:,2))];
znt_flx_new = znt(iflx_new,:);
znp_flx_new = znp(iflx_new,:);

% Construct new network: Z <= max_z and ( max(X) > min_xmf or flx_max > min_flx )
zn_flx_new = unique( [znt_flx_new;znp_flx_new], 'rows' );
zn_new = intersect(zn_zmax_new,union(zn_xmf_new,zn_flx_new,'rows'),'rows');

element = build_element_symbol;

% Build the species list (sunet)
nname_new = cell(size(zn_new,1),1);
for i = 1:size(zn_new,1)
    label = sprintf('%s%d',lower(strtrim(element(zn_new(i,1)+1,:))),zn_new(i,1)+zn_new(i,2));
    if strcmp(label,'n1')
        label = 'n';
    elseif strcmp(label,'h1')
        label = 'p';
    elseif strcmp(label,'h2')
        label = 'd';
    elseif strcmp(label,'h3')
        label = 't';
    end
    nname_new(i) = {label};
end

end

