function getnorm
global para
global norm
norm.a1=para.ang10r-para.ang10l;norm.a3=para.ang10peakmaxr-para.ang10peakmaxl;norm.a2=para.ang10l;norm.a4=para.ang10peakmaxl;
norm.b1=para.ang10raw;norm.b2= para.ang10mean;norm.b3=para.ang10max;
norm.c1=para.angpeak;norm.c2=para.angpeakloc;norm.c3= para.angpeaknum;
norm.d1=para.freq;norm.d2=para.freqpeaknum;norm.d3=para.freqmax;norm.d4=para.freqmaxloc;
norm.f1=para.tailpos;norm.f2= para.dur;norm.f3=para.marker;    norm.f4=para.length;
norm.g1=para.curvmax1;norm.g2=para.curvmaxloc1;norm.g3=para.curvmaxmax1;norm.g4=para.curvmaxmaxloc1;
norm.a1=mynorm(norm.a1);norm.a2=mynorm(norm.a2);norm.a3=mynorm(norm.a3);norm.a4=mynorm(norm.a4);
norm.b1=mynorm(norm.b1);norm.b2=mynorm(norm.b2);norm.b3=mynorm(norm.b3);
norm.c1=mynorm(norm.c1);norm.c2=mynorm(norm.c2);norm.c3=mynorm(norm.c3);
norm.d1=mynorm(norm.d1);norm.d2=mynorm(norm.d2);norm.d3=mynorm(norm.d3);norm.d4=mynorm(norm.d4);
norm.f1=mynorm(norm.f1);norm.f2=mynorm(norm.f2);
norm.g1=mynorm(norm.g1);norm.g2=mynorm(norm.g2);norm.g3=mynorm(norm.g3);norm.g4=mynorm(norm.g4);
end