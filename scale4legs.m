function [x,y] =  scale4legs( u, ur, v, vr)
%scale4legs   map one or two variables into range -1 to +1
% which is the range of the argument for legendre polynomials
%
%  x = scale4legs( u, ur)  maps the variable u which falls in the
%  range  ur = [ umin umax ] into the range [ -1 +1 ] using the
%  transformation
%                  x =  2*(u - umin)/(umax - umin)    - 1
%
%  Note: umax must be greater than umin.
%
%  [x,y] = scale4legs( u, ur, v, vr) performs scalings for both u
%  and v. This is a useful preliminary when constructing surfaces 
%  from a double legendre series.  vr = [vmin vmax]. v is mapped
%  onto y in the same way that u is mapped onto x.
%
%  see also unscale4legs


if length(ur)~=2 | ur(1)>=ur(2)
   error('ur has improper shape or contents')
end
i=find( u < ur(1) | u > ur(2));
if ~isempty(i)
   error(' u contains elements < umin or > umax')
end
x=u - ur(1);
x= 2*x./(ur(2)-ur(1))  -1;
if nargin==4
   if length(vr)~=2 | vr(1)>=vr(2)
      error('vr has improper shape or contents')
   end
   i=find( v < vr(1) | v > vr(2));
   if ~isempty(i)
      error(' v contains elements < vmin or > vmax')
   end
   y=v - vr(1);
   y= 2*y./(vr(2)-vr(1))  -1;
end

     