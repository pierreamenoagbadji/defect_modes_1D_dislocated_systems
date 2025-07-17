function val = cutoff(x, a, b, sigma, tol)
  % CUTOFF(x, a, b, sigma, tol) defines the usual cut-off function
  % which vanishes in x = a and x = b. sigma plays the role of variance.

  if (nargin < 5)
    tol = 1e-8;
  end

  if (nargin < 4)
    sigma = 1;
  end

  val = zeros(size(x));
  x0 = (x - a) / (b - a);
  I = (x0 > tol & x0 < 1 - tol);
  val(I) = exp((1.0 - 0.25 ./ (x0(I) .* (1 - x0(I)))) / (2*sigma^2));
end
