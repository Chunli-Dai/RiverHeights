function [p,pd] = legsd(n,x)
% LEGSD  LEGSD(n,x) evaluates legendre polynomials Pi(x) for i=0,1,..,n
% and their first derivatives with respect to x. 
% n is an integer, and x is a real scalar or vector. For each element
% of x, legsd computes the n+1 legendre polynomials of order 0 thru n 
% (returned in matrix p) and their 1st derivatives (returned in matrix pd). 
% If x has size [1,1], [1,m], or [m,1] then legsd returns
% three matrices of (same) size [1,n+1], [n+1,m], or [m,n+1] respectively.

% M.G. Bevis   10-9-90

if max(size(n)) ~= 1
   error('n must be scalar')
end
% need to check also n is positive integer, all elements of x are real.
% warn if x lies outside in range -1 to +1 

[a,b]= size(x);
if min(a,b) ~= 1
  error(' X must be a scalar or a vector')
end
xr=x(:)';         % xr is a row vector even if x is not
p=ones(n+1,max(a,b)); 
pd=zeros(n+1,max(a,b));
if n>0
   p(2,:)=xr;
   pd(2,:)=p(1,:);
   if n>1
      for l=2:n;
        i=l+1;   rl=1.0/l;
        p(i,:)= (2.0-rl)*xr.*p(i-1,:) - (1.0-rl)*p(i-2,:);
        pd(i,:)= l*p(i-1,:) + xr.*pd(i-1,:);
       end
   end
end
if b==1      % if x is a  column vector
   p=p'; pd=pd';
end