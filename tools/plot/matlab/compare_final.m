function [aerr,rerr] = compare_final(ts1,ts2,varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if ( nargin > 2 )
    verbose = varargin{1};
else
    verbose = false;
end

ny = ts1.ny;

enuc1 = cumsum(ts1.edot{:} .* ts1.tdel{:});
enuc2 = cumsum(ts2.edot{:} .* ts2.tdel{:});

nuc_name = build_isotope_symbol(ts1.zz{:},ts1.aa{:});

aerr = abs( ts2.xn{1}(:,end) - ts1.xn{1}(:,end) );
rerr = aerr ./ ts1.xn{1}(:,end);

[~,isort] = sort(rerr);
aerr(ny+1) = norm(aerr,2);
rerr(ny+1) = aerr(ny+1) ./ norm(ts1.xn{1}(:,end),2);

aerr(ny+2) = abs(ts2.t9{1}(end) - ts1.t9{1}(end));
rerr(ny+2) = aerr(ny+2) ./ ts1.t9{1}(end);

aerr(ny+3) = abs(enuc2(end) - enuc1(end));
rerr(ny+3) = aerr(ny+3) ./ enuc1(end);

if ( verbose )
    disp(sprintf("%5s %5s\t%10s\t%16s\t%16s\t%16s\t%16s",'i','isort','name','X1','X2','|dX|','|dX| / |X1|'));
    fmt = "%5d %5d\t%10s\t%16.8e\t%16.8e\t%16.8e\t%16.8e";
    for i = 1:ny
        disp(sprintf(fmt, i, isort(i), nuc_name{isort(i)}, ts1.xn{1}(isort(i),end), ts2.xn{1}(isort(i),end), aerr(isort(i)), rerr(isort(i))));
    end
    fmt = "%5s %5s\t%10s\t%16s\t%16s\t%16.8e\t%16.8e";
    disp(sprintf(fmt, '', '', '2-norm', '', '', aerr(ny+1), rerr(ny+1)));
    fmt = "%5s %5s\t%10s\t%16.8e\t%16.8e\t%16.8e\t%16.8e";
    disp(sprintf(fmt, '', '', 'T', ts1.t9{1}(end), ts2.t9{1}(end), aerr(ny+2), rerr(ny+2)));
    disp(sprintf(fmt, '', '', 'E_nuc', enuc1(end), enuc2(end), aerr(ny+3), rerr(ny+3)));
end

end

