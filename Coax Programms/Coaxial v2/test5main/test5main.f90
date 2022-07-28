    program Test5
    implicit complex*16(c,w,z), real*8(a-b,d-h,o-v,x-y)
    dimension ztemp(58), betam(14), qwork(2000), circle1(25), circle2(25)
    dimension z1(14), z2(14), z3(14), z(14), x2(14), x3(14), x4(14)
    dimension dC(4),dL(4), Um(2,2), aC(9), bC(4), dZ0(2,2), er(2)
    !--------------------------------------------------------------------------
    d1 = 2.; d2 = 5.; e = 10.;
    !--------------------------------------------------------------------------
    pi =  dacos(-1.d0)		  ! число пи
    data zi   /(0.d0, 1.d0) / ! мнимая единица
    !	Исходные данные:
    n = 1       ! количество линий
    nn = 14     ! Общее кол-во вершин в исх. многоугол. (вкл. доп.точки нуля и беcконечности)
    r1 = d1/2.
    r2 = d2/2.
    iprint = 0
    iguess = 1
    m=24;   ! m - кол-во углов по окружности сечения каждого провода СЛ.

    do i = 1, m+1
        alf = 2./m
        circle1(i) = r1*exp(-zi*(i+5)*alf*pi); ! Коорд.вершин n+1
    end do
    do i = 1, m+1
        circle2(i) = r2*exp(-zi*(i+5)*alf*pi); ! Коорд.вершин n+1
    end do 
   z1(1)  = circle1(13)
   z1(2)  = circle1(14) 
   z1(3)  = circle1(15) 
   z1(4)  = circle1(16)      
   z1(5)  = circle1(17)      
   z1(6)  = circle1(18)      
   z1(7)  = circle1(19)      
   z1(8)  = circle2(19)      
   z1(9)  = circle2(18)      
   z1(10) = circle2(17)      
   z1(11) = circle2(16)      
   z1(12) = circle2(15)      
   z1(13) = circle2(14)
   z1(14) = circle2(13)
   
    ! Углы при вершинах
    betam(1) =-0.5 + alf
    do i = 1, 5
        betam(i+1) = alf
    end do
    betam(7) =-0.5 + alf
    betam(8) =-0.5 - alf
    do i = 1, 5
        betam(i+8) = -alf
    end do
    betam(14) =-0.5 - alf
    ! Контрольная печать исходных данных
    do 11 i = 1,nn
11  write(6,100) i, z1(i)
100 format(' z1 =', i3, 2x, f6.3, 1x, f6.3)
    do 13 i = 1,nn
13  write(6,103) i, betam(i)
103 format(' betam =', i3, 2x, f6.3, 1x, f6.3)
    pause '----------------------------'

    z1c = dcmplx(r1, r1) ! Задание конформ. центра z10 в исх. многоугол. z1
    print *,z1c
    nptsq = 6
    call qinit(nn,betam,nptsq,qwork)

    do 1 k = 1,nn
        z(k) = exp(dcmplx(0.d0, k-nn))
        !		print*, k, '  ', z(k)
1   continue
    !	  pause '*****'

    tol = 1d-6	! 1.d-14

    call scsolv(iprint,iguess,tol,errest,nn,c,z,z1c,z1,betam,nptsq,qwork)
    !z1(k) - вершины многоугольника (в оригинале были обозначены как w(k) )
    !z(k) -  круговая каноническая область

    !Дробно-линейное отображение круга z(k) на верх.полуплоскость (окружность -> веществ.ось х2)
    z20 =zi;				     ! Задание конформ.центра z20 на верхней полуплоскости z2
    z201=-zi;						!z3
    do 2 k = 1, nn
        z2(k)=(z(k)*z201-z20)/(z(k)-1)  ! отображающая функция
        x2(k)=dreal(z2(k));			  ! выделение веществен.части (т.е. отбрасывание мнимой)
        print '(''z2 x2 '',i2, '' ('',f12.6,'','',f9.6'')'', f12.6)', k, z2(k), x2(k)
2   continue

    ! Нормировка к отрезку (-1...+1) на оси х3
    do 3 k = 1,nn
        x3(k) = 2*(x2(k) - x2(1)) / (x2(nn) - x2(1)) - 1.
        write(6,101) k, x3(k)
101     format(' x3 =', i3, 4x, f12.9 )
3   continue

    ! Выделение краевых точек электродов на оси х4
    x4(1)=x3(1);
    x4(2)=x3(7);
    x4(3)=x3(8);
    x4(4)=x3(14);
    print*,'----------------------------'
    !
    !	  do 5 i=1,2*n
    !	  write(6,102) i, x4(i)
    !5	  continue
    !102	  format(' x4 =', i3, 2x, f12.9)
    !	  print*,'----------------------------'
    !--------------------------------------------------------------------------
    !   n2 = 3
    !	M=1000
    !	call GHIONE(x4,aC,n2,M)

    ! Удаление "лишних" строк и столбцов, соответств. "лишним" электродам
    !	call refor(aC,bC,n)
    !	print*,'Capacitance matrix all in air bC'
    ! Матрица емкостей в воздухе
    !	call dprint(bC,n)
    !--------------------------------------------------------------------------
    ! Результирующая матрица емкостей (пФ/м)
    !	call capa(bC,dC,e,n)
    !   dC = e * bC;
    !	print*,'Capacitance matrix [C] (pF/m)'
    !	call dprint(dC,n)
    !--------------------------------------------------------------------------
    ! Результирующая матрица индуктивностей (нГн/м)
    !	call induc(bC,dL,n)
    !	print*,'Inductance matrix [L] (nH/m)'
    !	call dprint(dL,n)
    !--------------------------------------------------------------------------
    !    call Bet(dL, dC, dZ0, Um, er);
    !    print*,'Impedance matrix [Z] (Ohm)'
    !    call dprint(dZ0,n);
    !    print*,'Modal parameters [Um]'
    !    call dprint(Um,n);
    !    print*,'[er]'
    !    print*,er(1), er(2)
    pause
    end

