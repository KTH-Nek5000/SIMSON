% ***********************************************************************
%
% $HeadURL$
% $LastChangedDate$
% $LastChangedBy$
% $LastChangedRevision$
%
% ***********************************************************************
function [str2] = textline(text)
%        1234567890123456789012345678901234567890
  str = '                                        ';
  str2 = [str str];

  for i = 1:min(80,length(text))
    str2(i)=text(i);
  end