using Polynomials, LinearAlgebra

"""Implements the recurrence for the expected number of throws."""
function expected_throws(faces, levels)
  (e, d) = (ones(levels), zeros(levels))
  for k = 1:(levels - 1) # expectations
    for j = 1:min(k, faces)
        e[k + 1] = e[k + 1] + e[k - j + 1] / faces
    end
    e[k] = e[k]
  end
  for k = 1:levels # differences
    d[k] = e[k] - (2 * faces + 4) / (3 * (faces + 1))  - 2 * (k - 1) / (faces + 1)
  end
  return (e, d) # ingredients of table 1
end

"""Computes the second largest eigenvalue in magnitude."""
function second_eigenvalue(faces)
  n = faces
  coef = -ones(Int, n)
  push!(coef, n) 
  p = Polynomial(coef) / n # characteristic polynomial
  dp = derivative(p)
  root = roots(p) # eigenvalues
  d = norm.(root) # magnitude of eigenvalues
  perm = sortperm(d, rev = true) # sorted magnitudes
  z = root[perm[2]] # second largest eigenvalue
  z = conj(z) # complex conjugate
  y = exp(2 * pi * 1.0im / (n + 1)) # root of unity
  for k = 1:6 # Newton improvement
    y = y - p(y) / dp(y)
  end
  return (y, norm(p(y)), norm(z - y), n * norm(y) / (n - 1))
end

"""Main program."""
(faces, levels) = (6, 201);
(e, d) = expected_throws(faces, levels);
for level in [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 21, 51, 61, 81, 101, 201]
  println(level - 1," & ",round(e[level], sigdigits=5)," & ",
  round(d[level], sigdigits=5))
end
for faces in [2, 3, 4, 5, 6, 10, 20, 50, 100, 200]
  (u, a, b, c) = second_eigenvalue(faces)
  println(faces," & ",round(u, sigdigits=5)," & ",round(a, sigdigits=5)," & ",
  round(b, sigdigits=5)," & ",round(c, sigdigits=5))
end

