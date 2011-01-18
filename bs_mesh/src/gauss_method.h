#ifndef GAUSS_H
#define GAUSS_H

/*! \file gauss_method.cpp
	\brief solution linear system by gauss method
	\author Iskhakov Ruslan
	\date 2008-05-20 */

namespace grd_ecl
  {

  template <class fp_type_t>
  /*!	\brieve SWAP equations i && k
  	\param a[n*n] - matrix of left side
  	\param b[n*nn] - matrix of right side*/
  void reverse_eqns (fp_type_t *a, fp_type_t *b, const int n,const int nn, const int i,const  int k)
  {
    int j;
    fp_type_t tmp;

    for (j = 0; j < n; j++) //swap our string i&k
      {
        tmp = a[k*n + j];
        a[k*n + j] = a[i*n + j];
        a[i*n + j] = tmp;
      }

    for (j = 0; j < nn; j++)
      {
        tmp = b[k*nn+j];
        b[k*nn+j] = b[i*nn+j];
        b[i*nn+j] = tmp;
      }
  }

  /*! \brief solve linear system
  	\param a - matrix [may be modified in this function] n*n
  	\param b - output matrix [n*nn]*/
  template <class fp_type_t>
  int GaussSolve (fp_type_t *a /*matrix*/, fp_type_t *b, const int n, const int nn)
  {
    int i,j,k,k2,maxi;
    fp_type_t max;

    //first step of gauss
    for (k=0; k<n; k++)
      {
        max = fabs(a[k*n+k]);
        maxi = k;
        for (i=k+1; i<n; i++) //looking for maximum in column
          {
            if (fabs(a[i*n+k]) > max)
              {
                max = fabs(a[i*n+k]);
                maxi = i;
              }
          }
        if (max==0) return 1; //zero line

        if (maxi!=k) //if maximum on not his place
          reverse_eqns(a,b,n,nn,maxi,k);

        for (j=0; j<nn; j++)
          b[k*nn+j] /= a[k*n+k]; //right side vector

        for (j=k+1; j<n; j++)
          a[k*n+j] /= a[k*n+k]; // элементы с номерами столбцов меньшими к равны нулю

        a[k*n+k] = 1.0; // A[k][k]=1 нельзя засунуть в предыд. цикл т.к. после присвоения
        // нам понадобится его старое значение

        for (i=k+1; i<n; i++) //цикл по строкам
          {
            for (j=0;j<nn;j++)
              b[i*nn+j] -= b[k*nn+j]*a[i*n+k];

            for (j=k+1; j<n; j++) //вычитаем (к-тую) строчку, умноженную на ее элемент в к-том столбце
              a[i*n+j] -= a[k*n+j]*a[i*n+k];
            a[i*n+k] = 0.0;
          }
      }

    // ОБРАТНЫЙ ХОД МЕТОДА ГАУССА (явное решение системы)
    for (i = n - 2; i >= 0; i--)
      {
        for (j = i + 1; j < n; j++)
          for (k2 = 0; k2 < nn; k2++)
            b[i*nn+k2] -= a[i*n+j]*b[j*nn+k2];
      }
    return 0; //all it's ok
  }
}; //namespace grd_ecl
#endif //gauss_h
