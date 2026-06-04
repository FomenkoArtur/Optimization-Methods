function [x_m, f_m, t, n] = GoldenSection(a, b, E, f)
tic;
c = b - ((sqrt(5)-1)/2)*(b-a);
d = a + ((sqrt(5)-1)/2)*(b-a);
fc=f(c);
fd=f(d);
n=0;
while abs(b-a)>E
    n=n+1;
    if (fc<=fd)
        b=d;
        d=c;
        c=a+b-d;
        fd = fc;
        fc = f(c);
    else
        a=c;
        c=d;
        d=a+b-c;
        fc = fd;
        fd = f(d);
    end
end
x_m=(b+a)/2;
f_m=f(x_m);
t = toc;
end