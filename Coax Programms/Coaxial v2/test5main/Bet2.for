      subroutine Bet(dL, dC, dZ0, Um, er)
      real*8 bet2(2, 2), er(2), dC(2, 2), dL(2, 2), dZ0(2, 2), Um(2, 2)
      real*8 dlm(2,2), pi, d, dcs
      integer n
      dcs = 299792458;          ! Скорость света
      pi =  dacos(-1.d0)
      n = 2
      bet2 = matmul(dL, dC) * 1.e-21;
      er(1) = bet2(1,1) * dcs**2;
      er(2) = bet2(2,2) * dcs**2;
      bet2 = sqrt(bet2);
      call dminv(bet2, n, d);
      bet2(1,2) = 0;
      bet2(2,1) = 0;
      dC = dC * 1e-12;
      call EIGEN(bet2, Um, n, 0);
      dlm=matmul(matmul(dC, Um), bet2);
      call dminv(dlm, n, d);
      dZ0 = matmul(Um, dlm);
      return;
      end