function val = domainwall(x)

  % Attention: this way of writing the code
  % allows the user to evaluate the function
  % at x = ± inf

  Ipos = (x <= eps);
  Iint = (x >  eps & x <= 1-eps);
  
  val = zeros(size(x));

  val(Ipos) = 1;
  val(Iint) = (1 ./ (1 + exp((2*x(Iint) - 1)./(x(Iint).*(1 - x(Iint))))));
  val = 1 - val;

  % val = (x <= eps)...
  %     + (x >  eps & x <= 1-eps) .* (1 ./ (1 + exp((2*x - 1)./(x.*(1 - x))))) ;
end 