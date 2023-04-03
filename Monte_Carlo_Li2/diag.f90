subroutine diag(hm,n,e,v)
implicit none                 
real*8                                   :: tol
integer                                  :: n,ifail
real*8                                   :: hm(N,N),v(N,N),e(N)
real*8,dimension(:),allocatable          :: bet

allocate(bet(N))
  tol = 2.0d0**(-50)
  call f01ajf(n, tol, hm, n, e, bet, v, n)
  tol = 2.0d0**(-55)
  ifail = 1
  call f02amf(n, tol, e, bet, v, n, ifail)
deallocate(bet)
end subroutine diag

subroutine diag1(hm,om,n,e,v)
  implicit none                 
  integer                                  :: n,ifail
  real*8                                   :: hm(N,N),v(N,N),e(N),om(n,n)
  real*8,dimension(:),allocatable          :: el,dl

  allocate(el(n),dl(n))

  call f02aef(hm,n,om,n,n,e,v,n,dl,el,ifail)

  deallocate(el,dl)

end subroutine diag1

      subroutine f01aef(n,a,ia,b,ib,dl,ifail)
!     mark 2 release. nag copyright 1972
!     mark 3 revised.
!     mark 4.5 revised

!     reduc1
!     reduction of the general symmetric eigenvalue problem a * x =
!     lambda * b * x, with symmetric matrix a and symmetric
!     positive definite matrix b, to the equivalent standard
!     problem p * z = lambda * z. the upper triangle, including
!     diagonal elements, of a and b are given in the arrays a(n,n)
!     and b(n,n). l(b = l * lt) is formed in the
!     remaining strictly lower triangle of the array b with its
!     diagonal elements in the array dl(n) , and the lower
!     triangle of the symmetric matrix p(p = inv(l) * a * inv( lt))
!     is formed in the lower triangle of the array a, including the
!     diagonal elements. hence the diagonal elements of a are lost.
!     the subroutine will fail if b, perhaps on account of rounding
!     errors, is not positive definite - sets ifail = 0 if
!     successful else ifail = 1.
!     1st december 1971

      integer isave, ifail, i, n, i1, j, kk, k, j1, ia, ib, p01aaf
      double precision srname
      double precision x, y, a(ia,n), b(ib,n), dl(n)
      data srname /8h f01aef /
      isave = ifail
      do 100 i=1,n
         i1 = i - 1
         do 80 j=i,n
            x = b(i,j)
            if (i1.eq.0) go to 40
            do 20 kk=1,i1
               k = i1 - kk + 1
               x = x - b(i,k)*b(j,k)
   20       continue
   40       if (i.ne.j) go to 60
            if (x.le.0.0d0) go to 320
            y = dsqrt(x)
            dl(i) = y
            go to 80
   60       b(j,i) = x/y
   80    continue
  100 continue
!     l has been formed in array b
      do 180 i=1,n
         y = dl(i)
         i1 = i - 1
         do 160 j=i,n
            x = a(i,j)
            if (i1.eq.0) go to 140
            do 120 kk=1,i1
               k = i1 - kk + 1
               x = x - b(i,k)*a(j,k)
  120       continue
  140       a(j,i) = x/y
  160    continue
  180 continue
!     the transpose of the upper triangle of
!     inv(l) * a has been formed in the lower
!     triangle of array a
      do 300 j=1,n
         j1 = j - 1
         do 280 i=j,n
            x = a(i,j)
            i1 = i - 1
            if (i1.lt.j) go to 220
            do 200 kk=j,i1
               k = i1 - kk + j
               x = x - a(k,j)*b(i,k)
  200       continue
  220       if (j1.eq.0) go to 260
            do 240 kk=1,j1
               k = j1 - kk + 1
               x = x - a(j,k)*b(i,k)
  240       continue
  260       a(i,j) = x/dl(i)
  280    continue
  300 continue
      ifail = 0
      return
  320 ifail = p01aaf(isave,1,srname)
      return
      end
      subroutine f01aff(n, im1, im2, b, ib, dl, z, iz)
!     mark 2 release. nag copyright 1972
!     mark 4.5 revised

!     rebaka
!     this subroutine performs, on the matrix of eigenvectors, z,
!     stored
!     in columns im1 to im2 of the array z(n,n), a backward
!     substitution
!     lt* x = z, over- writing x on z. the diagonal elements of l
!     must be stored in the array dl(n), and the remaining
!     triangle in the strictly lower triangle of the array b(n,n).
!     the subroutines f01aef and f01bdf leave l in this
!     desired form. if x denotes any column of the resultant matrix
!     x, then x satisfies xt * b * x = zt * z, where b = l * lt.
!     1st august 1971

      integer j, im1, im2, ii, n, i, i2, k, ib, iz
      double precision x, b(ib,n), dl(n), z(iz,im2)
      do 80 j=im1,im2
         do 60 ii=1,n
            i = n - ii + 1
            x = z(i,j)
            i2 = i + 1
            if (i2.gt.n) go to 40
            do 20 k=i2,n
               x = x - b(k,i)*z(k,j)
   20       continue
   40       z(i,j) = x/dl(i)
   60    continue
   80 continue
      return
      end
      subroutine f01ajf(n, atol, a, ia, d, e, z, iz)
!     mark 2 release. nag copyright 1972
!     mark 4 revised.
!     mark 4.5 revised
!     mark 5c revised

!     tred2
!     this subroutine reduces the given lower triangle of a
!     symmetric matrix, a, stored in the array a(n,n), to
!     tridiagonal form using householders reduction. the diagonal
!     of the result is stored in the array d(n) and the
!     sub-diagonal in the last n - 1 stores of the array e(n)
!     (with the additional element e(1) = 0). the transformation
!     matrices are accumulated in the array z(n,n). the array
!     a is left unaltered unless the actual parameters
!     corresponding to a and z are identical.
!     1st august 1971

      integer i, ia, ii, iz, j1, j, k, l, n
      double precision atol, f, g, h, hh, a(ia,n), d(n), e(n), z(iz,n)
      do 40 i=1,n
         do 20 j=1,i
            z(i,j) = a(i,j)
   20    continue
   40 continue
      if (n.eq.1) go to 280
      do 260 ii=2,n
         i = n - ii + 2
         l = i - 2
         f = z(i,i-1)
         g = 0.0d0
         if (l.eq.0) go to 80
         do 60 k=1,l
            g = g + z(i,k)*z(i,k)
   60    continue
   80    h = g + f*f
!     if g is too small for orthogonality to be
!     guaranteed the transformation is skipped
         if (g.gt.atol) go to 100
         e(i) = f
         h = 0.0d0
         go to 240
  100    l = l + 1
         g = dsqrt(h)
         if (f.ge.0.0d0) g = -g
         e(i) = g
         h = h - f*g
         z(i,i-1) = f - g
         f = 0.0d0
         do 180 j=1,l
            z(j,i) = z(i,j)/h
            g = 0.0d0
!     form element of a*u
            do 120 k=1,j
               g = g + z(j,k)*z(i,k)
  120       continue
            j1 = j + 1
            if (j1.gt.l) go to 160
            do 140 k=j1,l
               g = g + z(k,j)*z(i,k)
  140       continue
!     form element of p
  160       e(j) = g/h
            f = f + g*z(j,i)
  180    continue
!     form k
         hh = f/(h+h)
!     form reduced a
         do 220 j=1,l
            f = z(i,j)
            g = e(j) - hh*f
            e(j) = g
            do 200 k=1,j
               z(j,k) = z(j,k) - f*e(k) - g*z(i,k)
  200       continue
  220    continue
  240    d(i) = h
  260 continue
  280 e(1) = 0.0d0
      d(1) = 0.0d0
!     accumulation of transformation matrices
      do 400 i=1,n
         l = i - 1
         if (d(i).eq.0.0d0) go to 360
         do 340 j=1,l
            g = 0.0d0
            do 300 k=1,l
               g = g + z(i,k)*z(k,j)
  300       continue
            do 320 k=1,l
               z(k,j) = z(k,j) - g*z(k,i)
  320       continue
  340    continue
  360    d(i) = z(i,i)
         z(i,i) = 1.0d0
         if (l.eq.0) go to 400
         do 380 j=1,l
            z(i,j) = 0.0d0
            z(j,i) = 0.0d0
  380    continue
  400 continue
      return
      end

      subroutine f02aef(a,ia,b,ib,n,r,v,iv,dl,e,ifail)
!     mark 2 release. nag copyright 1972
!     mark 3 revised.
!     mark 4.5 revised

!     eigenvalues and eigenvectors of a-lambda*b
!     1st december 1971

      integer p01aaf,isave,ifail,n,ia,ib,iv
      double precision srname
      double precision tol,a(ia,n),b(ib,n),r(n),v(iv,n),dl(n),e(n)
      data srname /8h f02aef /
      isave = ifail
      ifail = 1
      call f01aef(n, a, ia, b, ib, dl, ifail)
      if (ifail.eq.0) go to 20
      ifail = p01aaf(isave,ifail,srname)
      return
   20 tol = 2.0d0**(-50)
      call f01ajf(n, tol, a, ia, r, e, v, iv)
      tol = 2.0d0**(-55)
      ifail = 1
      call f02amf(n, tol, r, e, v, iv, ifail)
      if (ifail.eq.0) go to 40
      ifail = p01aaf(isave,2,srname)
      return
   40 call f01aff(n, 1, n, b, ib, dl, v, iv)
      return
      end
      subroutine f02amf(n,acheps,d,e,z,iz,ifail)
!     mark 2 release. nag copyright 1972
!     mark 3 revised.
!     mark 4 revised.
!     mark 4.5 revised

!     tql2
!     this subroutine finds the eigenvalues and eigenvectors of a
!     tridiagonal matrix, t, given with its diagonal elements in
!     the array d(n) and its sub-diagonal elements in the last n
!     - 1 stores of the array e(n), using ql transformations. the
!     eigenvalues are overwritten on the diagonal elements in the
!     array d in ascending order. the eigenvectors are formed in
!     the array z(n,n), overwriting the accumulated
!     transformations as supplied by the subroutine f01ajf. the
!     subroutine will fail if any one eigenvalue takes more than 30
!     iterations.
!     1st april 1972
!
      integer p01aaf, isave, ifail, n, i, l, j, m, i1, m1, ii, k, iz
      double precision srname
      double precision b, f, h, acheps, g, p, r, c, s, d(n), e(n), z(iz,n)
      data srname /8h f02amf /
      isave = ifail
      if (n.eq.1) go to 40
      do 20 i=2,n
         e(i-1) = e(i)
   20 continue
   40 e(n) = 0.0d0
      b = 0.0d0
      f = 0.0d0
      do 300 l=1,n
         j = 0
         h = acheps*(dabs(d(l))+dabs(e(l)))
         if (b.lt.h) b = h
!     look for small sub-diag element
         do 60 m=l,n
            if (dabs(e(m)).le.b) go to 80
   60    continue
   80    if (m.eq.l) go to 280
  100    if (j.eq.30) go to 400
         j = j + 1
!     form shift
         g = d(l)
         h = d(l+1) - g
         if (dabs(h).ge.dabs(e(l))) go to 120
         p = h*0.5d0/e(l)
         r = dsqrt(p*p+1.0d0)
         h = p + r
         if (p.lt.0.0d0) h = p - r
         d(l) = e(l)/h
         go to 140
  120    p = 2.0d0*e(l)/h
         r = dsqrt(p*p+1.0d0)
         d(l) = e(l)*p/(1.0d0+r)
  140    h = g - d(l)
         i1 = l + 1
         if (i1.gt.n) go to 180
         do 160 i=i1,n
            d(i) = d(i) - h
  160    continue
  180    f = f + h
!     ql transformation
         p = d(m)
         c = 1.0d0
         s = 0.0d0
         m1 = m - 1
         do 260 ii=l,m1
            i = m1 - ii + l
            g = c*e(i)
            h = c*p
            if (dabs(p).lt.dabs(e(i))) go to 200
            c = e(i)/p
            r = dsqrt(c*c+1.0d0)
            e(i+1) = s*p*r
            s = c/r
            c = 1.0d0/r
            go to 220
  200       c = p/e(i)
            r = dsqrt(c*c+1.0d0)
            e(i+1) = s*e(i)*r
            s = 1.0d0/r
            c = c/r
  220       p = c*d(i) - s*g
            d(i+1) = h + s*(c*g+s*d(i))
!     form vector
            do 240 k=1,n
               h = z(k,i+1)
               z(k,i+1) = s*z(k,i) + c*h
               z(k,i) = c*z(k,i) - s*h
  240       continue
  260    continue
         e(l) = s*p
         d(l) = c*p
         if (dabs(e(l)).gt.b) go to 100
  280    d(l) = d(l) + f
  300 continue
!     order eigenvalues and eigenvectors
      do 380 i=1,n
         k = i
         p = d(i)
         i1 = i + 1
         if (i1.gt.n) go to 340
         do 320 j=i1,n
            if (d(j).ge.p) go to 320
            k = j
            p = d(j)
  320    continue
  340    if (k.eq.i) go to 380
         d(k) = d(i)
         d(i) = p
         do 360 j=1,n
            p = z(j,i)
            z(j,i) = z(j,k)
            z(j,k) = p
  360    continue
  380 continue
      ifail = 0
      return
  400 ifail = p01aaf(isave,1,srname)
      return
      end

      

      subroutine p01aaz
!     routine to stop program
      write(6,*)'stop'
      end




      integer function p01aaf(ifail, error, srname)
!     mark 1 release.  nag copyright 1971
!     mark 3 revised
!     mark 4a revised, ier-45
!     mark 4.5 revised
!     mark 7 revised (dec 1978)
!     returns the value of error or terminates the program.
      integer error, ifail, nout

      double precision srname
!*** for v5 compatability only
!     test if no error detected
      if (error.eq.0) go to 20
      nout = 9
!     test for soft failure
      if (mod(ifail,10).eq.1) go to 10
!     hard failure
      write (nout,99999) srname, error
!     stopping mechanism may also differ
      call p01aaz
!     stop
!     soft fail
!     test if error messages suppressed
   10 if (mod(ifail/10,10).eq.0) go to 20
      write (nout,99999) srname, error
   20 p01aaf = error
      return
99999 format (1h0, 38herror detected by nag library routine , a8,11h - ifail = , i5//)
      end


