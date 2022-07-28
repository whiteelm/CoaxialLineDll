        SUBROUTINE DPRINT(A,N)
C  ÐÀÑÏÅ×ÀÒÊÀ ÂÅÙÅÑÒÂÅÍÍÛÕ ÊÂÀÄÐÀÒÍÛÂÕ ÌÀÒÐÈÖ ÐÀÇÌÅÐÎÌ ÍÅ ÁÎËÅÅ 24x24
        real*8 A(1)
        IF(N.GT.24)print*, '*** DPRINT N>24'
        IF(N.GT.24)RETURN
        JN1=8
        JN2=16
        JN3=24
        IF(N.LE.24) KK=3
        IF(N.LE.24) JN3=N
        IF(N.LE.16)  KK=2
        IF(N.LE.16)  JN2=N
        IF(N.LE.8)  KK=1
        IF(N.LE.8)  JN1=N
        DO1  K=1,KK
        WRITE(*,'(/)')
        DO1  I=1,N
        IF(K.EQ.1)WRITE(*,'(1x,8f10.4)')(A((j-1)*n+i),J=1,JN1)
        IF(K.EQ.2)WRITE(*,'(1x,8f10.4)')(A((j-1)*n+i),J=9,JN2)
        IF(K.EQ.3)WRITE(*,'(1x,8f10.4)')(A((j-1)*n+i),J=17,JN3)
1       CONTINUE
        WRITE(*,'(/)')
        RETURN
        END