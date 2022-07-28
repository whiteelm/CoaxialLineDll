    program CoaxialConsole
    real*8 pi, dL, dC, dZ0, a, b, er
    integer mod
    pi=3.14159265359
    dL=0
    dC=0
    dZ0=138.0922
    b=10
    a=1
    er=1
    mod=3
    if (mod == 1) then
        dL=2*10e-7*log(b/a);
        dC=(2*pi*8.85*10e-12*er)/log(b/a);
        dZ0=sqrt(dL/dC);
    else if (mod == 2) then   
        b = a/(10**(-(dZ0*sqrt(er))/138));
        dL=2*10e-7*log(b/a);
        dC=(2*pi*8.85*10e-12*er)/log(b/a);
    else if (mod == 3) then
        a = b*10**(-(dZ0*sqrt(er))/138);
        dL=2*10e-7*log(b/a);
        dC=(2*pi*8.85*10e-12*er)/log(b/a);
    endif
    print *, 'mod', mod
    print *, 'dL', dL*1e6
    print *, 'dC', dC*1e12
    print *, 'dZ0', dZ0
    print *, 'a', a
    print *, 'b', b
    print *, 'er', er
    pause
    end program CoaxialConsole

