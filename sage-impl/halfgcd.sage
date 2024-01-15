from sage.all import *

# taken from: https://scrapbox.io/crypto-writeup-public/half-GCD

MIN = 10
def pdivmod(u, v):
  q = u // v
  r = u - q*v
  return (q, r)

def hgcd(u, v, min_degree=MIN):
  x = u.parent().gen()

  if u.degree() < v.degree():
    u, v = v, u

  if 2*v.degree() < u.degree() or u.degree() < min_degree:
    q = u // v
    return matrix([[1, -q], [0, 1]])

  m = u.degree() // 2
  b0 = u // (x^m)
  b1 = v // (x^m)

  R = hgcd(b0, b1, min_degree=min_degree)
  DE = R * matrix([[u], [v]])
  d, e = DE[0,0], DE[1,0]
  q, f = pdivmod(d, e)

  g0 = e // x^(m//2)
  g1 = f // x^(m//2)

  S = hgcd(g0, g1, min_degree=min_degree)
  return S * matrix([[0, 1], [1, -q]]) * R

def poly_gcd_euclid(u: Polynomial, v: Polynomial):
  while v != 0:
    u, v = v, u % v
  return u.monic()

def poly_gcd(u: Polynomial, v: Polynomial):
  if u.degree() < v.degree():
    u, v = v, u

  while True:
    if v == 0:
      return u.monic()

    if u % v == 0:
      return v.monic()

    if u.degree() < MIN:
      return poly_gcd_euclid(u, v)

    R = hgcd(u, v)
    B = R * matrix([[u], [v]])
    b0, b1 = B[0,0], B[1,0]
    u, v = b1, b0 % b1
