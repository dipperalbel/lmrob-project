c Goto's BLAS at least needs a XERBLA
      subroutine xerbla(srname, info)
      character*6 srname
      integer info
      end

      subroutine test1(iflag)
      double complex zx(2), ztemp, zres, zdotu
      integer iflag
      zx(1) = (3.1d0,1.7d0)
      zx(2) = (1.6d0,-0.6d0)
      zres = zdotu(2, zx, 1, zx, 1)
      ztemp = (0.0d0,0.0d0)
      do 10 i = 1,2
 10      ztemp = ztemp + zx(i)*zx(i)
      if(abs(zres - ztemp) > 1.0d-10) then
        iflag = 1
      else
        iflag = 0
      endif
      end

