function h = fspecials(type,P1,P2,varargin)


switch length(P1)
case 1
h=fspecial1(type,P1,P2(1,1)); % needs scalar std.deviation
case 2
h=fspecial2s(type,P1,P2,varargin{:});
case 3
h=fspecial3s(type,P1,P2);
end;
