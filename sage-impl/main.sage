from time import time
from typing import List
from sage.all import *

load("halfgcd.sage")

p = 2**31-1
G = Zmod(p)
F.<x> = PolynomialRing(G)

def reduce(polys: List[Polynomial], mat: Matrix):
  """
  不明な項を 1 次増やす代わりに次数を k-1 次減らす
  計算量: O(k^2d)
  """
  k = len(polys)

  start_deg = polys[0].degree()
  for i in range(k):
    if polys[i].degree() < start_deg:
      polys[i] += polys[0]
  assert all([poly.degree() == start_deg for poly in polys])

  # O(d)
  def reduce_poly(i, j):
    q = polys[i] // polys[j]
    assert q.degree() <= 1
    mat[i] -= q * mat[j] # O(k*d/k)
    polys[i] -= q * polys[j] # O(d)

  for divisor_ind in range(k-1):
    for dividend_ind in range(k):
      if divisor_ind == dividend_ind: continue
      reduce_poly(dividend_ind, divisor_ind)

  reduce_poly(-2, -1)
  mat[-2] += mat[-1]
  polys[-2] += polys[-1]

  end_deg = polys[0].degree()
  assert start_deg - end_deg == k - 1
  assert all([poly.degree() == end_deg for poly in polys])

def calc_M(polys: List):
  """
  M * vector(polys) の要素の次数が ceil(d / k) 次程度になるような行列 M を返す
  計算量: O(k^2 dlog^2d)
  """
  d, k = polys[0].degree(), len(polys)

  if polys[0].degree() <= k + 1:
    return matrix.identity(F, k)

  M_1 = calc_M([poly // x^(d // 2) for poly in polys])

  new_polys = list((M_1 * matrix(polys).T).T[0])  # O(k^2 dlogd)

  if new_polys[0].degree() <= k + 1:
    return M_1
  
  reduce(new_polys, M_1) # O(k^2d)

  if new_polys[0].degree() <= k + 1:
    return M_1

  reduced_deg = d - new_polys[0].degree()
  gained_pollution = ceil(reduced_deg / (k - 1))

  M_2 = calc_M([poly // x^gained_pollution for poly in new_polys])
  return M_2 * M_1 # O(k^3 d/klog(d/k))

def apply_M(polys: List, M: Matrix):
  """
  M * vector(polys) を計算し、そのうち二項を返す
  計算量: O(k(d/k+k)log(d/k+k))
  """
  required_degree = max([max([poly.degree() for poly in row]) for row in M]) + len(polys)
  
  M2 = M[:2]
  fs = [poly[:required_degree+1] for poly in polys]
  V2 = M2 * matrix(fs).T

  return V2[0, 0][:required_degree+1], V2[1, 0][:required_degree+1]

d = 1000
k = 10

fs = [F.random_element(d) for _ in range(k)]
a = G.random_element()
cs = [f(a) for f in fs]

vs = [f - c for f, c in zip(fs, cs)]

start = time()
half_gcd_res = poly_gcd(vs[0], vs[1])
half_gcd_elapsed = time() - start
print(f'[*] {half_gcd_elapsed=:.5}')

start = time()
M = calc_M(vs)
calc_M_elapsed = time() - start
print(f'[*] {calc_M_elapsed=:.5}')

start = time()
u, v = apply_M(vs, M)
apply_M_elapsed = time() - start
print(f'[*] {apply_M_elapsed=:.5}')

start = time()
res = poly_gcd(u, v)
compute_gcd_elapsed = time() - start
print(f'[*] {compute_gcd_elapsed=:.5}')

print()
print(f'[+] {x-a=}')
print(f'[+] {half_gcd_res=}')
print(f'[+] {res=}')

assert x-a == res == half_gcd_res
