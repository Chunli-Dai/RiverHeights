function p = legs(n,x)
% LEGS  LEGS(n,x) evaluates legendre polynomials Pi(x) for i=0,1,..,n. 
% n is an integer, and x is a real scalar or vector. For each element
% of x, legs computes the n+1 legendre polynomials of order 0 thru n.
% If x has size [1,1], [1,m], or [m,1] then legs returns
% a real matrix of size [1,n+1], [n+1,m], or [m,n+1] respectively.

% M.G. Bevis   5-9-90

if max(size(n)) ~= 1
  error('n must be scalar')
end
% need to check also n is positive integer, all elements of x are real.
% warn if x lies outside in range -1 to +1 

[a,b]= size(x);
if min(a,b) ~= 1
  error(' X must be a scalar or a vector')
end
xr=x(:)';        % xr is a row vector even if x is a column vector
p=ones(n+1,max(a,b));    % two rows needed for storing various orders of Pn(x)
if n>0
  p(2,:)=xr;
  if n>1
    for l=2:n;
      i=l+1; rl=1.0/l;
      p(i,:)= (2.0-rl)*xr.*p(i-1,:) - (1.0-rl)*p(i-2,:);
    end     
  end
end
if b==1           % x is a column vector
  p=p';    
end         