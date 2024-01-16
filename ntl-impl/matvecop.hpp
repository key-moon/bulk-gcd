#include <NTL/lzz_pX.h>
#include <NTL/vector.h>
#include <NTL/matrix.h>
#include <NTL/BasicThreadPool.h>

NTL_CLIENT

const int MMMUL_THRESHOLD = 10000;
const int MVMUL_THRESHOLD = 10000;
const int VVADD_THRESHOLD = 10000;

template <typename T=zz_pX, typename TContext=zz_pContext, bool Parallel=false>
class MatVecOp {
  public:
  static void mul(Vec<T>& c, const Mat<T>& a, const Vec<T>& b) {
    assert(a.NumCols() == b.length());
    int n = a.NumRows(), m = a.NumCols();

    c.SetLength(n);

    #define LOOP_1(first, last) \
      for (int i = first; i < last; i++) { \
        c[i] = 0; \
        for (int j = 0; j < m; j++) { \
          c[i] += a[i][j] * b[j];\
        }\
      }
    if (Parallel and MVMUL_THRESHOLD < deg(a[0][0])) {
      TContext context; context.save();
      NTL_EXEC_RANGE(n, first, last)
      context.restore();
      LOOP_1(first, last)
      NTL_EXEC_RANGE_END
    }
    else {
      LOOP_1(0, n)
    }
  }

  static void mul(Mat<T>& c, const Mat<T>& a, const Mat<T>& b) {
    assert(a.NumCols() == b.NumRows());
    int n = a.NumRows(), m = a.NumCols(), l = b.NumCols();

    c.SetDims(n, l);

    #define LOOP_2(first, last) \
      for (int i = first; i < last; i++) { \
        for (int j = 0; j < l; j++) { \
          c[i][j] = 0; \
          for (int k = 0; k < m; k++) { \
            c[i][j] += a[i][k] * b[k][j]; \
          } \
        } \
      }

    if (Parallel and MMMUL_THRESHOLD < deg(a[0][0]) * n) {
      TContext context; context.save();
      NTL_EXEC_RANGE(n, first, last)
      context.restore();
      LOOP_2(first, last)
      NTL_EXEC_RANGE_END
    }
    else {
      LOOP_2(0, n)
    }
  }

  template<typename U>
  static void mul_add_inplace(Vec<T>& a, const Vec<T>& b, U q) {
    assert(a.length() == b.length());
    int n = a.length();

    #define LOOP_3(first, last) \
      for (int i = first; i < last; i++) { \
        a[i] += b[i] * q; \
      }
    if (Parallel and VVADD_THRESHOLD < deg(b[0])) {
      TContext context; context.save();
      NTL_EXEC_RANGE(n, first, last)
      context.restore();
      LOOP_3(first, last)
      NTL_EXEC_RANGE_END
    }
    else {
      LOOP_3(0, n)
    }
  }
};
