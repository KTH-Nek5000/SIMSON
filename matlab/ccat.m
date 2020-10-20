% ***********************************************************************
%
% $HeadURL$
% $LastChangedDate$
% $LastChangedBy$
% $LastChangedRevision$
%
% ***********************************************************************
function res=ccat(dim,varargin);
res=[];
for ind=1:length(varargin);
  res=cat(dim,res,varargin{ind});
end
