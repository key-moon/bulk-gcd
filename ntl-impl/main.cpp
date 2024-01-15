#include <bits/stdc++.h>
#include <NTL/vector.h>
#include <NTL/matrix.h>
#include <NTL/lzz_pX.h>

#include "matvecop.hpp"

#ifndef LOG_NAME
#define LOG_NAME "log.csv"
#endif

#ifdef DO_PARALLEL
using Op=MatVecOp<zz_pX, zz_pContext, true>;
#else
using Op=MatVecOp<zz_pX, zz_pContext, false>;
#endif

NTL_CLIENT

void reduce(Vec<zz_pX>& polys, Mat<zz_pX>& m) {
  int k = polys.length();
  int start_deg = deg(polys[0]);
  for (int i = 0; i < k; i++) {
    if (deg(polys[i]) < start_deg) {
      polys[i] += polys[0];
    }
  }

  auto reduce_poly = [&](int i, int j) {
    zz_pX q = polys[i] / polys[j];
    assert(deg(q) <= 1);
    Op::mul_add_inplace(m[i], m[j], -q);
    polys[i] -= q * polys[j];
  };

  for (int divisor_ind = 0; divisor_ind < k - 1; divisor_ind++) {
    for (int dividend_ind = 0; dividend_ind < k; dividend_ind++) {
      if (divisor_ind == dividend_ind) continue;
      reduce_poly(dividend_ind, divisor_ind);
    }
  }

  reduce_poly(k - 2, k - 1);
  Op::mul_add_inplace(m[k - 2], m[k - 1], zz_p(1));
  polys[k - 2] += polys[k - 1];
  
  int end_deg = deg(polys[0]);
  assert(start_deg - end_deg == k - 1);
}

void _bulk_gcd(Mat<zz_pX>& m, Vec<zz_pX>& polys) {
  int k = polys.length();
  int d = deg(polys[0]);

  if (d <= k + 1) {
    m.SetDims(k, k);
    for (int i = 0; i < k; i++) {
      m[i][i] = 1;
    }
    return;
  }

  int mid = d / 2;
  Mat<zz_pX> m1;
  {
    Vec<zz_pX> arg_polys;
    arg_polys.SetLength(k);
    for (int i = 0; i < k; i++) {
      arg_polys[i] = polys[i] >> mid;
    }
    _bulk_gcd(m1, arg_polys);
  }

  Vec<zz_pX> new_polys;
  {
    int required_degree = d - (mid / k) * (k - 1) + k + 1;

    if (required_degree <= k + 1) { m = m1; return; }

    for (auto&& poly : polys) {
       if (required_degree < deg(poly)) poly.SetLength(required_degree + 1);
    }
    
    Op::mul(new_polys, m1, polys);

    for (auto&& poly : new_polys) {
      if (required_degree < deg(poly)) {
        while (poly[required_degree] == 0) {
          poly.SetLength(deg(poly));
          required_degree--;
        }
        poly.SetLength(required_degree + 1);
      }
    }

    if (deg(new_polys[0]) <= k + 1) { m = m1; return; }
  }

  reduce(new_polys, m1);

  if (deg(new_polys[0]) <= k + 1) { m = m1; return; }

  Mat<zz_pX> m2;
  {
    int reduced_deg = d - deg(new_polys[0]);
    int gained_pollution = (reduced_deg + (k - 1) - 1) / (k - 1); // ceil(reduced_deg / (k - 1))
    Vec<zz_pX> arg_polys;
    arg_polys.SetLength(k);
    for (int i = 0; i < k; i++) {
      arg_polys[i] = new_polys[i] >> gained_pollution;
    }
    _bulk_gcd(m2, arg_polys);
  }

  Op::mul(m, m2, m1);
}

void bulk_gcd(Mat<zz_pX>& m, const Vec<zz_pX>& polys) {
  Vec<zz_pX> _polys = polys;
  _bulk_gcd(m, _polys);
}

void apply_m(Vec<zz_pX>& res, Vec<zz_pX>& polys, Mat<zz_pX>& m) {
  int max_deg = 0;
  for (int i = 0; i < m.NumRows(); i++) {
    for (int j = 0; j < m.NumCols(); j++) {
      int deg = ::deg(m.get(i, j));
      if (max_deg < deg) max_deg = deg;
    }
  }

  int required_degree = max_deg + polys.length() + 10;

  Mat<zz_pX> m2;
  m2.SetDims(2, m.NumCols());
  m2[0] = m[0]; m2[1] = m[1];

  for (auto&& f : polys) {
    if (required_degree < deg(f)) {
      f.SetLength(required_degree + 1);
    }
  }

  Op::mul(res, m2, polys);

  if (required_degree < deg(res[0])) {
    while (res[0][required_degree] == 0) required_degree--;
    res[0].SetLength(required_degree + 1);
    res[1].SetLength(required_degree + 1);
  }
}

int main() {
  setbuf(stdout, NULL);
  SetNumThreads(8);

  zz_p::init((1LL << 31) - 1);
  zz_p a;
  random(a);

  double start;
  vector<int> k_cands { 3, 5, 10, 15, 20 };

  int iter = 5;
  for (int n = 5000; n <= 200000; n += 5000) {
    for (int kind = 0; kind < k_cands.size(); kind++) {
      int k = k_cands[kind];
      vector<vector<double>> res_bucket(5);
      for (int cnt = 0; cnt < iter; cnt++) {
        Vec<zz_pX> polys;
        polys.SetLength(k);
        for (int i = 0; i < k; i++) {
          random(polys[i], n);
          polys[i] -= eval(polys[i], a);
        }

        zz_pX res_poly_gcd;
        start = GetTime();
        GCD(res_poly_gcd, polys[0], polys[1]);
        res_bucket[0].emplace_back(GetTime() - start);

        Mat<zz_pX> m;
        start = GetTime();
        bulk_gcd(m, polys);
        res_bucket[1].emplace_back(GetTime() - start);

        Vec<zz_pX> reduced;
        start = GetTime();
        apply_m(reduced, polys, m);
        res_bucket[2].emplace_back(GetTime() - start);

        zz_pX res_poly;
        start = GetTime();
        GCD(res_poly, reduced[0], reduced[1]);
        res_bucket[3].emplace_back(GetTime() - start);
      }

      vector<double> res(5);
      for (int i = 0; i < 4; i++) {
        sort(res_bucket[i].begin(), res_bucket[i].end());
        for (int j = 1; j < res_bucket[i].size() - 1; j++) {
          res[i] += res_bucket[i][j];
        }
        res[i] /= res_bucket[i].size() - 2;
      }
      res[4] += res[2] + res[3];

      FILE *fp;
      fp = fopen(LOG_NAME, "a");
      fprintf(fp, "%d,%d", n, k);
      for (int i = 0; i < res.size(); i++) {
        fprintf(fp, ",%f", res[i]);
      }
      fprintf(fp, "\n");
      fclose(fp);
    }
  }
}
