    subroutine Coaxial(er, a, b, dC, dL, dZ0, mod)
    !dec$ attributes dllexport ::COAXIAL
    !DEC$ ATTRIBUTES VALUE :: mod
    !DEC$ ATTRIBUTES REFERENCE :: er, a, b, dC, dL, em, dZ0
    real*8 pi, dL, dC, dZ0, a, b, er
    integer mod;
    pi=3.14159265359;
    if (mod == 1) then
        dL=.2*log(b/a);
        dC=2*pi*8.85*er/log(b/a);
        dZ0=138*log10(b/a)/sqrt(er);
    else 
        a = b*10**(-(dZ0*sqrt(er))/138);
        dL=.2*10e-7*log(b/a);
        dC=.2*pi*8.85*er/log(b/a);
    endif
    end