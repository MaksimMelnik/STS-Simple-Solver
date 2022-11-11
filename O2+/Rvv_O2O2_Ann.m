function out=Rvv_O2O2_Ann(n, T, O2)
% R VV with k_i->i-1 ^j-1->j and reverse function like in Annušová 2018 
% paper
% 23.02.2021
k = 1.3807e-23;
ABdata1=[
1 2 2.2802E-17 1.2679E+00 12 38 7.5015E-30 5.5905E+00
1 3 2.2105E-17 1.3096E+00 12 39 1.0320E-30 5.8269E+00
1 4 1.6591E-17 1.3860E+00 12 40 1.1033E-31 6.0965E+00
1 5 1.0611E-17 1.4740E+00 12 41 1.3597E-32 6.3371E+00
1 6 7.3215E-18 1.5446E+00 13 14 1.1715E-14 1.0580E+00
1 7 3.6462E-18 1.6540E+00 13 15 7.9854E-15 1.1207E+00
1 8 1.8430E-18 1.7567E+00 13 16 3.7123E-15 1.2372E+00
1 9 8.7833E-19 1.8674E+00 13 17 1.9428E-15 1.3284E+00
1 10 4.0885E-19 1.9786E+00 13 18 9.3864E-16 1.4315E+00
1 11 1.7580E-19 2.0985E+00 13 19 4.2131E-16 1.5410E+00
1 12 6.4292E-20 2.2438E+00 13 20 1.9069E-16 1.6481E+00
1 13 2.8884E-20 2.3541E+00 13 21 9.1580E-17 1.7464E+00
1 14 1.3949E-20 2.4534E+00 13 22 3.7319E-17 1.8673E+00
1 15 5.5817E-21 2.5812E+00 13 23 1.3727E-17 1.9956E+00
1 16 3.3050E-21 2.6426E+00 13 24 4.5044E-18 2.1476E+00
1 17 1.4746E-21 2.7501E+00 13 25 1.3833E-18 2.3010E+00
1 18 5.1214E-22 2.8930E+00 13 26 3.9193E-19 2.4663E+00
1 19 1.5975E-22 3.0535E+00 13 27 1.1819E-19 2.6217E+00
1 20 9.8884E-23 3.1047E+00 13 28 2.6524E-20 2.8137E+00
1 21 3.7198E-23 3.2308E+00 13 29 6.3801E-21 2.9907E+00
1 22 1.4741E-23 3.3507E+00 13 30 1.1025E-21 3.2195E+00
1 23 4.5105E-24 3.5051E+00 13 31 2.5499E-22 3.4025E+00
1 24 1.6581E-24 3.6378E+00 13 32 4.9321E-23 3.6119E+00
1 25 5.5900E-25 3.7753E+00 13 33 7.7191E-24 3.8445E+00
1 26 1.9113E-25 3.9107E+00 13 34 1.2130E-24 4.0669E+00
1 27 6.6239E-26 4.0441E+00 13 35 1.1560E-25 4.3573E+00
1 28 1.7785E-26 4.2082E+00 13 36 5.4020E-27 4.7639E+00
1 29 6.5429E-27 4.3203E+00 13 37 4.0695E-28 5.0950E+00
1 30 1.3880E-27 4.5191E+00 13 38 1.6398E-29 5.5145E+00
1 31 2.1911E-28 4.7584E+00 13 39 1.5843E-30 5.8029E+00
1 32 6.2638E-29 4.9097E+00 13 40 1.2888E-31 6.1112E+00
1 33 1.0904E-29 5.1241E+00 13 41 1.0678E-32 6.4102E+00
1 34 1.8886E-30 5.3317E+00 14 15 2.3455E-14 9.8648E-01
1 35 1.9558E-31 5.6218E+00 14 16 2.1343E-14 1.0035E+00
1 36 1.9120E-32 5.9157E+00 14 17 1.2097E-14 1.0904E+00
1 37 1.4568E-33 6.2485E+00 14 18 6.0355E-15 1.1908E+00
1 38 8.6605E-35 6.6136E+00 14 19 2.9213E-15 1.2944E+00
1 39 7.6659E-36 6.9118E+00 14 20 1.2194E-15 1.4166E+00
1 40 3.6090E-37 7.3013E+00 14 21 5.5779E-16 1.5255E+00
1 41 7.8562E-39 7.8033E+00 14 22 2.3216E-16 1.6447E+00
2 3 1.2973E-16 1.1702E+00 14 23 6.9082E-17 1.8065E+00
2 4 1.0618E-16 1.2361E+00 14 24 2.8113E-17 1.9267E+00
2 5 6.6533E-17 1.3310E+00 14 25 8.2721E-18 2.0879E+00
2 6 2.9463E-17 1.4684E+00 14 26 2.4435E-18 2.2476E+00
2 7 1.8861E-17 1.5403E+00 14 27 7.2721E-19 2.4056E+00
2 8 1.0295E-17 1.6321E+00 14 28 1.5997E-19 2.6012E+00
2 9 5.6629E-18 1.7199E+00 14 29 4.2116E-20 2.7652E+00
2 10 2.7628E-18 1.8274E+00 14 30 6.7800E-21 3.0033E+00
2 11 1.6269E-18 1.9008E+00 14 31 1.5377E-21 3.1883E+00
2 12 4.9839E-19 2.0756E+00 14 32 2.6278E-22 3.4206E+00
2 13 2.2672E-19 2.1834E+00 14 33 4.9540E-23 3.6206E+00
2 14 9.1913E-20 2.3078E+00 14 34 3.1566E-24 3.9800E+00
2 15 4.7712E-20 2.3944E+00 14 35 3.4415E-25 4.2563E+00
2 16 1.6169E-20 2.5442E+00 14 36 2.6872E-26 4.5769E+00
2 17 7.2929E-21 2.6499E+00 14 37 1.0567E-27 5.0099E+00
2 18 3.2226E-21 2.7587E+00 14 38 4.2051E-29 5.4273E+00
2 19 1.0983E-21 2.9041E+00 14 39 4.3601E-30 5.6975E+00
2 20 4.8590E-22 3.0078E+00 14 40 3.0662E-31 6.0306E+00
2 21 1.6446E-22 3.1532E+00 14 41 1.7828E-32 6.3755E+00
2 22 7.5853E-23 3.2510E+00 15 16 2.2813E-14 1.0100E+00
2 23 2.3610E-23 3.4056E+00 15 17 1.8398E-14 1.0459E+00
2 24 1.0079E-23 3.5119E+00 15 18 1.0367E-14 1.1302E+00
2 25 3.7281E-24 3.6337E+00 15 19 6.6741E-15 1.1915E+00
2 26 1.2558E-24 3.7714E+00 15 20 2.5502E-15 1.3306E+00
2 27 4.2912E-25 3.9068E+00 15 21 1.2076E-15 1.4337E+00
2 28 1.1321E-25 4.0735E+00 15 22 5.1125E-16 1.5511E+00
2 29 3.3419E-26 4.2196E+00 15 23 1.5399E-16 1.7118E+00
2 30 7.6177E-27 4.4049E+00 15 24 6.3033E-17 1.8316E+00
2 31 1.5060E-27 4.6139E+00 15 25 1.8589E-17 1.9929E+00
2 32 1.6713E-28 4.9037E+00 15 26 5.4848E-18 2.1532E+00
2 33 3.0314E-29 5.1143E+00 15 27 1.6279E-18 2.3119E+00
2 34 5.9092E-30 5.3070E+00 15 28 3.9462E-19 2.4933E+00
2 35 7.4109E-31 5.5640E+00 15 29 8.5131E-20 2.6885E+00
2 36 1.0798E-31 5.7973E+00 15 30 1.3954E-20 2.9266E+00
2 37 9.9597E-33 6.1009E+00 15 31 3.0850E-21 3.1186E+00
2 38 4.0843E-34 6.5212E+00 15 32 5.5943E-22 3.3356E+00
2 39 3.3504E-35 6.8354E+00 15 33 7.8417E-23 3.5849E+00
2 40 3.7969E-36 7.0981E+00 15 34 6.8395E-24 3.8958E+00
2 41 1.7051E-37 7.4874E+00 15 35 5.5497E-25 4.2186E+00
3 4 6.6026E-16 1.0511E+00 15 36 4.5816E-26 4.5314E+00
3 5 6.5477E-16 1.0838E+00 15 37 3.5659E-27 4.8557E+00
3 6 3.2187E-16 1.2131E+00 15 38 1.0409E-28 5.3193E+00
3 7 2.1474E-16 1.2825E+00 15 39 6.3642E-30 5.6709E+00
3 8 1.4179E-16 1.3504E+00 15 40 6.6079E-31 5.9439E+00
3 9 6.4350E-17 1.4735E+00 15 41 5.3942E-32 6.2376E+00
3 10 3.2148E-17 1.5766E+00 16 17 3.7777E-14 9.6012E-01
3 11 1.2409E-17 1.7185E+00 16 18 2.7181E-14 1.0128E+00
3 12 6.1585E-18 1.8204E+00 16 19 2.0924E-14 1.0506E+00
3 13 3.1478E-18 1.9123E+00 16 20 1.4221E-14 1.1043E+00
3 14 1.3455E-18 2.0322E+00 16 21 5.8144E-15 1.2342E+00
3 15 5.1176E-19 2.1673E+00 16 22 2.1534E-15 1.3736E+00
3 16 3.1051E-19 2.2281E+00 16 23 8.9405E-16 1.4885E+00
3 17 1.2670E-19 2.3479E+00 16 24 3.0770E-16 1.6354E+00
3 18 5.5529E-20 2.4588E+00 16 25 8.3849E-17 1.8095E+00
3 19 2.0800E-20 2.5911E+00 16 26 2.4809E-17 1.9702E+00
3 20 9.6515E-21 2.6919E+00 16 27 7.3389E-18 2.1301E+00
3 21 3.2255E-21 2.8372E+00 16 28 1.5867E-18 2.3294E+00
3 22 1.5881E-21 2.9228E+00 16 29 2.9236E-19 2.5507E+00
3 23 6.6049E-22 3.0293E+00 16 30 4.6649E-20 2.7912E+00
3 24 1.4676E-22 3.2378E+00 16 31 1.2527E-20 2.9537E+00
3 25 4.7412E-23 3.3818E+00 16 32 2.3082E-21 3.1695E+00
3 26 1.5551E-23 3.5234E+00 16 33 3.1444E-22 3.4233E+00
3 27 5.1768E-24 3.6628E+00 16 34 2.4303E-23 3.7543E+00
3 28 1.3234E-24 3.8340E+00 16 35 3.4242E-24 3.9882E+00
3 29 2.5853E-25 4.0404E+00 16 36 1.8013E-25 4.3714E+00
3 30 4.6065E-26 4.2630E+00 16 37 1.2500E-26 4.7141E+00
3 31 9.4196E-27 4.4602E+00 16 38 2.5293E-28 5.2319E+00
3 32 1.6205E-27 4.6839E+00 16 39 1.9537E-29 5.5432E+00
3 33 2.5463E-28 4.9129E+00 16 40 1.4655E-30 5.8590E+00
3 34 3.1646E-29 5.1698E+00 16 41 8.9984E-32 6.1930E+00
3 35 3.1242E-30 5.4691E+00 17 18 4.2360E-14 9.6084E-01
3 36 4.4194E-31 5.7090E+00 17 19 2.6422E-14 1.0347E+00
3 37 5.1128E-32 5.9781E+00 17 20 1.6148E-14 1.1053E+00
3 38 4.7097E-33 6.2771E+00 17 21 9.8116E-15 1.1760E+00
3 39 3.2316E-34 6.6129E+00 17 22 4.5475E-15 1.2829E+00
3 40 2.1223E-35 6.9610E+00 17 23 1.6869E-15 1.4165E+00
3 41 3.0889E-36 7.1747E+00 17 24 6.2307E-16 1.5526E+00
4 5 1.1364E-15 1.0492E+00 17 25 1.8809E-16 1.7123E+00
4 6 8.0223E-16 1.1244E+00 17 26 5.5977E-17 1.8728E+00
4 7 3.9643E-16 1.2419E+00 17 27 1.6577E-17 2.0329E+00
4 8 2.6312E-16 1.3109E+00 17 28 3.9361E-18 2.2185E+00
4 9 1.8223E-16 1.3676E+00 17 29 6.6937E-19 2.4523E+00
4 10 8.1117E-17 1.4965E+00 17 30 1.4446E-19 2.6478E+00
4 11 4.1528E-17 1.5938E+00 17 31 3.4606E-20 2.8248E+00
4 12 2.0290E-17 1.7007E+00 17 32 5.0106E-21 3.0798E+00
4 13 8.4029E-18 1.8248E+00 17 33 6.9531E-22 3.3318E+00
4 14 3.4316E-18 1.9500E+00 17 34 8.0661E-23 3.6029E+00
4 15 1.3092E-18 2.0847E+00 17 35 7.1999E-24 3.9055E+00
4 16 6.2920E-19 2.1815E+00 17 36 3.6418E-25 4.2929E+00
4 17 2.2163E-19 2.3246E+00 17 37 2.4856E-26 4.6380E+00
4 18 1.1967E-19 2.4022E+00 17 38 5.3289E-28 5.1466E+00
4 19 4.9551E-20 2.5196E+00 17 39 3.2688E-29 5.4966E+00
4 20 1.6921E-20 2.6638E+00 17 40 2.6565E-30 5.7960E+00
4 21 6.7079E-21 2.7856E+00 17 41 1.4974E-31 6.1493E+00
4 22 2.8100E-21 2.8994E+00 18 19 5.4054E-14 9.4350E-01
4 23 8.7712E-22 3.0491E+00 18 20 4.2884E-14 9.7871E-01
4 24 3.7445E-22 3.1595E+00 18 21 2.5021E-14 1.0594E+00
4 25 1.1996E-22 3.3048E+00 18 22 1.2645E-14 1.1557E+00
4 26 3.9028E-23 3.4477E+00 18 23 4.7049E-15 1.2913E+00
4 27 1.2890E-23 3.5882E+00 18 24 2.3639E-15 1.3829E+00
4 28 3.2093E-24 3.7641E+00 18 25 5.9394E-16 1.5707E+00
4 29 6.0619E-25 3.9769E+00 18 26 1.7919E-16 1.7301E+00
4 30 1.2254E-25 4.1854E+00 18 27 5.3398E-17 1.8901E+00
4 31 3.9271E-26 4.3141E+00 18 28 1.1506E-17 2.0910E+00
4 32 4.8284E-27 4.5930E+00 18 29 3.3471E-18 2.2444E+00
4 33 8.4685E-28 4.8082E+00 18 30 5.9870E-19 2.4691E+00
4 34 7.7024E-29 5.1150E+00 18 31 8.1495E-20 2.7361E+00
4 35 1.0384E-29 5.3610E+00 18 32 1.4512E-20 2.9541E+00
4 36 1.4332E-30 5.6028E+00 18 33 1.9670E-21 3.2114E+00
4 37 1.2670E-31 5.9062E+00 18 34 1.3648E-22 3.5585E+00
4 38 9.8740E-33 6.2303E+00 18 35 1.1899E-23 3.8722E+00
4 39 1.3162E-33 6.4693E+00 18 36 8.8208E-25 4.1983E+00
4 40 1.0756E-34 6.7812E+00 18 37 6.0093E-26 4.5413E+00
4 41 7.9532E-36 7.1029E+00 18 38 1.3134E-27 5.0464E+00
5 6 1.5087E-15 1.0706E+00 18 39 6.7888E-29 5.4209E+00
5 7 8.9477E-16 1.1642E+00 18 40 5.2658E-30 5.7314E+00
5 8 7.2814E-16 1.2026E+00 18 41 4.6030E-31 6.0113E+00
5 9 4.0578E-16 1.2978E+00 19 20 6.9309E-14 9.2426E-01
5 10 1.9357E-16 1.4135E+00 19 21 5.5716E-14 9.5945E-01
5 11 1.1260E-16 1.4919E+00 19 22 3.2310E-14 1.0382E+00
5 12 5.0052E-17 1.6143E+00 19 23 1.4496E-14 1.1474E+00
5 13 2.0922E-17 1.7377E+00 19 24 5.9195E-15 1.2721E+00
5 14 6.2641E-18 1.9118E+00 19 25 1.8495E-15 1.4295E+00
5 15 2.9666E-18 2.0133E+00 19 26 5.7122E-16 1.5865E+00
5 16 1.3448E-18 2.1181E+00 19 27 1.7241E-16 1.7455E+00
5 17 4.5131E-19 2.2681E+00 19 28 3.7330E-17 1.9465E+00
5 18 2.5508E-19 2.3366E+00 19 29 9.9189E-18 2.1135E+00
5 19 1.0527E-19 2.4547E+00 19 30 1.4310E-18 2.3697E+00
5 20 3.7953E-20 2.5927E+00 19 31 3.0205E-19 2.5674E+00
5 21 1.6123E-20 2.7057E+00 19 32 4.1522E-20 2.8285E+00
5 22 6.8612E-21 2.8177E+00 19 33 5.8500E-21 3.0799E+00
5 23 2.5599E-21 2.9408E+00 19 34 6.9241E-22 3.3419E+00
5 24 9.0153E-22 3.0801E+00 19 35 5.6250E-23 3.6675E+00
5 25 2.8646E-22 3.2267E+00 19 36 2.8325E-24 4.0561E+00
5 26 9.2468E-23 3.3708E+00 19 37 1.7720E-25 4.4132E+00
5 27 3.0296E-23 3.5126E+00 19 38 4.5835E-27 4.8932E+00
5 28 7.5919E-24 3.6869E+00 19 39 3.5580E-28 5.2021E+00
5 29 1.7382E-24 3.8708E+00 19 40 1.8255E-29 5.5747E+00
5 30 4.8218E-25 4.0284E+00 19 41 9.4460E-31 5.9381E+00
5 31 8.2356E-26 4.2592E+00 20 21 9.9021E-14 8.8850E-01
5 32 2.2307E-26 4.4162E+00 20 22 7.0640E-14 9.4063E-01
5 33 3.7810E-27 4.6370E+00 20 23 3.5503E-14 1.0390E+00
5 34 5.8058E-28 4.8621E+00 20 24 1.6370E-14 1.1449E+00
5 35 5.0581E-29 5.1766E+00 20 25 5.6200E-15 1.2901E+00
5 36 6.7508E-30 5.4219E+00 20 26 1.8006E-15 1.4430E+00
5 37 6.7689E-31 5.7102E+00 20 27 5.5520E-16 1.5998E+00
5 38 6.8504E-32 5.9926E+00 20 28 1.2172E-16 1.8000E+00
5 39 6.9838E-33 6.2731E+00 20 29 3.4908E-17 1.9540E+00
5 40 5.0287E-34 6.6087E+00 20 30 6.9113E-18 2.1675E+00
5 41 4.4733E-35 6.9038E+00 20 31 7.6360E-19 2.4677E+00
6 7 1.3607E-15 1.1305E+00 20 32 1.5252E-19 2.6699E+00
6 8 8.9603E-16 1.2045E+00 20 33 1.8599E-20 2.9405E+00
6 9 5.3127E-16 1.2909E+00 20 34 2.5599E-21 3.1872E+00
6 10 3.2073E-16 1.3715E+00 20 35 1.7162E-22 3.5343E+00
6 11 1.3114E-16 1.5050E+00 20 36 8.1649E-24 3.9317E+00
6 12 8.4843E-17 1.5697E+00 20 37 5.6048E-25 4.2757E+00
6 13 3.5720E-17 1.6924E+00 20 38 1.6985E-26 4.7303E+00
6 14 9.9154E-18 1.8769E+00 20 39 7.1149E-28 5.1382E+00
6 15 5.3064E-18 1.9620E+00 20 40 2.7023E-29 5.5506E+00
6 16 2.2450E-18 2.0758E+00 20 41 2.1206E-30 5.8473E+00
6 17 1.0731E-18 2.1737E+00 21 22 1.1119E-13 8.8472E-01
6 18 4.1557E-19 2.3021E+00 21 23 6.1509E-14 9.6949E-01
6 19 1.6138E-19 2.4285E+00 21 24 3.2777E-14 1.0579E+00
6 20 6.8961E-20 2.5399E+00 21 25 1.1946E-14 1.1957E+00
6 21 2.8135E-20 2.6589E+00 21 26 3.9662E-15 1.3444E+00
6 22 1.1696E-20 2.7740E+00 21 27 1.2498E-15 1.4988E+00
6 23 3.6082E-21 2.9257E+00 21 28 3.0055E-16 1.6859E+00
6 24 1.5267E-21 3.0375E+00 21 29 5.9468E-17 1.8974E+00
6 25 4.8321E-22 3.1848E+00 21 30 1.1213E-17 2.1177E+00
6 26 1.5533E-22 3.3295E+00 21 31 2.2952E-18 2.3244E+00
6 27 5.0685E-23 3.4720E+00 21 32 3.4542E-19 2.5707E+00
6 28 1.2631E-23 3.6471E+00 21 33 4.2994E-20 2.8397E+00
6 29 3.5689E-24 3.8000E+00 21 34 5.8396E-21 3.0889E+00
6 30 1.1851E-24 3.9299E+00 21 35 4.7409E-22 3.4043E+00
6 31 2.3632E-25 4.1379E+00 21 36 1.6875E-23 3.8470E+00
6 32 3.6285E-26 4.3799E+00 21 37 1.0644E-24 4.2040E+00
6 33 6.1119E-27 4.6016E+00 21 38 3.3183E-26 4.6540E+00
6 34 1.1952E-27 4.7945E+00 21 39 1.3828E-27 5.0631E+00
6 35 1.2915E-28 5.0704E+00 21 40 6.2316E-29 5.4543E+00
6 36 1.4834E-29 5.3407E+00 21 41 4.1579E-30 5.7752E+00
6 37 1.9472E-30 5.5903E+00 22 23 1.2941E-13 8.7505E-01
6 38 1.4402E-31 5.9235E+00 22 24 8.1119E-14 9.4233E-01
6 39 1.4565E-32 6.2024E+00 22 25 3.3557E-14 1.0642E+00
6 40 1.5555E-33 6.4769E+00 22 26 1.1986E-14 1.2040E+00
6 41 1.2531E-34 6.7908E+00 22 27 3.9464E-15 1.3532E+00
7 8 2.8725E-15 1.0754E+00 22 28 9.0499E-16 1.5488E+00
7 9 2.2208E-15 1.1272E+00 22 29 1.7871E-16 1.7633E+00
7 10 1.7325E-15 1.1765E+00 22 30 4.0659E-17 1.9525E+00
7 11 9.5321E-16 1.2689E+00 22 31 5.2082E-18 2.2297E+00
7 12 5.8322E-16 1.3463E+00 22 32 1.1306E-18 2.4242E+00
7 13 3.0215E-16 1.4410E+00 22 33 1.3830E-19 2.6963E+00
7 14 1.5721E-16 1.5339E+00 22 34 1.0229E-20 3.0380E+00
7 15 6.4695E-17 1.6581E+00 22 35 9.8904E-22 3.3292E+00
7 16 1.9975E-17 1.8255E+00 22 36 5.0707E-23 3.7154E+00
7 17 9.7553E-18 1.9211E+00 22 37 2.7985E-24 4.0921E+00
7 18 3.8843E-18 2.0467E+00 22 38 5.5647E-26 4.6091E+00
7 19 1.4270E-18 2.1839E+00 22 39 4.0917E-27 4.9262E+00
7 20 4.6955E-19 2.3342E+00 22 40 2.6036E-28 5.2664E+00
7 21 3.1522E-19 2.3759E+00 22 41 1.3210E-29 5.6299E+00
7 22 1.0054E-19 2.5313E+00 23 24 2.0239E-13 8.2477E-01
7 23 2.7121E-20 2.7055E+00 23 25 1.1208E-13 9.0870E-01
7 24 1.2616E-20 2.8015E+00 23 26 4.6896E-14 1.0285E+00
7 25 3.8997E-21 2.9526E+00 23 27 1.6938E-14 1.1664E+00
7 26 1.2235E-21 3.1012E+00 23 28 3.9147E-15 1.3641E+00
7 27 3.8982E-22 3.2473E+00 23 29 8.9346E-16 1.5614E+00
7 28 9.4328E-23 3.4270E+00 23 30 1.4408E-16 1.8068E+00
7 29 1.7605E-23 3.6381E+00 23 31 4.0381E-17 1.9625E+00
7 30 4.6284E-24 3.8020E+00 23 32 5.4320E-18 2.2275E+00
7 31 9.4821E-25 4.0027E+00 23 33 6.5250E-19 2.5032E+00
7 32 2.7024E-25 4.1569E+00 23 34 7.2670E-20 2.7815E+00
7 33 5.2957E-26 4.3472E+00 23 35 4.8049E-21 3.1383E+00
7 34 4.9331E-27 4.6520E+00 23 36 2.2002E-22 3.5370E+00
7 35 5.3384E-28 4.9294E+00 23 37 1.2580E-23 3.9080E+00
7 36 4.1627E-29 5.2609E+00 23 38 2.2105E-25 4.4428E+00
7 37 3.4655E-30 5.5782E+00 23 39 1.3051E-26 4.7950E+00
7 38 2.8325E-31 5.8922E+00 23 40 1.1109E-27 5.0864E+00
7 39 3.4807E-32 6.1409E+00 23 41 3.8241E-29 5.5007E+00
7 40 5.0883E-33 6.3698E+00 24 25 2.3821E-13 8.1152E-01
7 41 9.1978E-34 6.5540E+00 24 26 1.2183E-13 9.0522E-01
8 9 5.1224E-15 1.0346E+00 24 27 4.9210E-14 1.0289E+00
8 10 5.1589E-15 1.0480E+00 24 28 1.6610E-14 1.1709E+00
8 11 3.4259E-15 1.1174E+00 24 29 2.5862E-15 1.4267E+00
8 12 2.3991E-15 1.1771E+00 24 30 4.9122E-16 1.6480E+00
8 13 1.3088E-15 1.2660E+00 24 31 1.1465E-16 1.8367E+00
8 14 5.8312E-16 1.3815E+00 24 32 1.8453E-17 2.0730E+00
8 15 2.6869E-16 1.4937E+00 24 33 2.1939E-18 2.3511E+00
8 16 1.2296E-16 1.6016E+00 24 34 2.5451E-19 2.6275E+00
8 17 5.7352E-17 1.7065E+00 24 35 2.3273E-20 2.9333E+00
8 18 2.1500E-17 1.8411E+00 24 36 6.7129E-22 3.4006E+00
8 19 7.8001E-18 1.9800E+00 24 37 3.8137E-23 3.7732E+00
8 20 2.8134E-18 2.1181E+00 24 38 6.5492E-25 4.3109E+00
8 21 1.0236E-18 2.2575E+00 24 39 5.2269E-26 4.6219E+00
8 22 4.5100E-19 2.3631E+00 24 40 2.5761E-27 4.9959E+00
8 23 1.5832E-19 2.4983E+00 24 41 1.2901E-28 5.3637E+00
8 24 5.5280E-20 2.6375E+00 25 26 3.3754E-13 7.7229E-01
8 25 1.6829E-20 2.7912E+00 25 27 1.7485E-13 8.6360E-01
8 26 5.1994E-21 2.9423E+00 25 28 6.8170E-14 9.8806E-01
8 27 1.6302E-21 3.0910E+00 25 29 1.9251E-14 1.1553E+00
8 28 3.8650E-22 3.2739E+00 25 30 3.8320E-15 1.3691E+00
8 29 9.5452E-23 3.4450E+00 25 31 7.1918E-16 1.5918E+00
8 30 1.7467E-23 3.6636E+00 25 32 9.2485E-17 1.8669E+00
8 31 5.5635E-24 3.7984E+00 25 33 1.0935E-17 2.1469E+00
8 32 1.0302E-24 4.0156E+00 25 34 5.7239E-19 2.5434E+00
8 33 2.1010E-25 4.2041E+00 25 35 3.6935E-20 2.9032E+00
8 34 1.6536E-26 4.5337E+00 25 36 3.0553E-21 3.2126E+00
8 35 2.0208E-27 4.7903E+00 25 37 1.6736E-22 3.5909E+00
8 36 1.4098E-28 5.1357E+00 25 38 2.7810E-24 4.1329E+00
8 37 1.1448E-29 5.4564E+00 25 39 1.9813E-25 4.4563E+00
8 38 5.7625E-31 5.8419E+00 25 40 5.5101E-27 4.9204E+00
8 39 5.6249E-32 6.1308E+00 25 41 3.9065E-28 5.2336E+00
8 40 1.0026E-32 6.3212E+00 26 27 4.8856E-13 7.2863E-01
8 41 8.7511E-34 6.6203E+00 26 28 2.0431E-13 8.4807E-01
9 10 7.2131E-15 1.0193E+00 26 29 7.3942E-14 9.7899E-01
9 11 5.8342E-15 1.0625E+00 26 30 1.4782E-14 1.1970E+00
9 12 4.4488E-15 1.1112E+00 26 31 2.2413E-15 1.4558E+00
9 13 2.4142E-15 1.2044E+00 26 32 4.7036E-16 1.6565E+00
9 14 1.2869E-15 1.2946E+00 26 33 5.5899E-17 1.9371E+00
9 15 6.9753E-16 1.3814E+00 26 34 3.7873E-18 2.2938E+00
9 16 2.7632E-16 1.5136E+00 26 35 2.4234E-19 2.6516E+00
9 17 1.5595E-16 1.5893E+00 26 36 1.4477E-20 3.0173E+00
9 18 7.1193E-17 1.6939E+00 26 37 7.6496E-22 3.4014E+00
9 19 2.8775E-17 1.8188E+00 26 38 1.4300E-23 3.9252E+00
9 20 8.7223E-18 1.9848E+00 26 39 9.7134E-25 4.2574E+00
9 21 3.5242E-18 2.1064E+00 26 40 2.7490E-26 4.7119E+00
9 22 9.6878E-19 2.2848E+00 26 41 1.2951E-27 5.0897E+00
9 23 2.9044E-19 2.4427E+00 27 28 6.2206E-13 7.0082E-01
9 24 1.2505E-19 2.5508E+00 27 29 2.1469E-13 8.4292E-01
9 25 3.7808E-20 2.7057E+00 27 30 6.2913E-14 1.0048E+00
9 26 1.1593E-20 2.8581E+00 27 31 1.3822E-14 1.2066E+00
9 27 3.6072E-21 3.0081E+00 27 32 2.4039E-15 1.4432E+00
9 28 8.4634E-22 3.1927E+00 27 33 2.9160E-16 1.7222E+00
9 29 1.7964E-22 3.3882E+00 27 34 3.2687E-17 2.0058E+00
9 30 5.6993E-23 3.5278E+00 27 35 1.7940E-18 2.3879E+00
9 31 8.2050E-24 3.7839E+00 27 36 7.1416E-20 2.8146E+00
9 32 1.8013E-24 3.9734E+00 27 37 3.6433E-21 3.2045E+00
9 33 2.7307E-25 4.2106E+00 27 38 6.2843E-23 3.7399E+00
9 34 2.5343E-26 4.5090E+00 27 39 2.8978E-24 4.1296E+00
9 35 3.5761E-27 4.7464E+00 27 40 1.3092E-25 4.5143E+00
9 36 2.7040E-28 5.0754E+00 27 41 5.9915E-27 4.8857E+00
9 37 2.1543E-29 5.3990E+00 28 29 8.4393E-13 6.6268E-01
9 38 1.0025E-30 5.7961E+00 28 30 2.4361E-13 8.3374E-01
9 39 8.2370E-32 6.1128E+00 28 31 7.4800E-14 9.8934E-01
9 40 1.2707E-32 6.3248E+00 28 32 1.7650E-14 1.1795E+00
9 41 1.1566E-33 6.6186E+00 28 33 2.2618E-15 1.4525E+00
10 11 7.9500E-15 1.0386E+00 28 34 1.3752E-16 1.8304E+00
10 12 5.6057E-15 1.0976E+00 28 35 9.7223E-18 2.1795E+00
10 13 3.5585E-15 1.1668E+00 28 36 5.1582E-19 2.5616E+00
10 14 2.0114E-15 1.2499E+00 28 37 2.6373E-20 2.9514E+00
10 15 9.6719E-16 1.3543E+00 28 38 3.9111E-22 3.5092E+00
10 16 4.0142E-16 1.4799E+00 28 39 2.7262E-23 3.8334E+00
10 17 1.5681E-16 1.6123E+00 28 40 1.2407E-24 4.2140E+00
10 18 6.6972E-17 1.7292E+00 28 41 2.4189E-26 4.7167E+00
10 19 2.4229E-17 1.8688E+00 29 30 1.3232E-12 6.0407E-01
10 20 9.9672E-18 1.9863E+00 29 31 6.0384E-13 7.0660E-01
10 21 4.1031E-18 2.1065E+00 29 32 1.7417E-13 8.7276E-01
10 22 1.6626E-18 2.2267E+00 29 33 2.5407E-14 1.1295E+00
10 23 4.9149E-19 2.3862E+00 29 34 1.7686E-15 1.4884E+00
10 24 2.0122E-19 2.5037E+00 29 35 1.2447E-16 1.8388E+00
10 25 6.0635E-20 2.6592E+00 29 36 6.3784E-18 2.2314E+00
10 26 1.8529E-20 2.8122E+00 29 37 2.9992E-19 2.6359E+00
10 27 5.7425E-21 2.9629E+00 29 38 4.2113E-21 3.2020E+00
10 28 1.3403E-21 3.1483E+00 29 39 5.4458E-22 3.4269E+00
10 29 3.3616E-22 3.3198E+00 29 40 1.4374E-23 3.8903E+00
10 30 7.5726E-23 3.5118E+00 29 41 5.9139E-25 4.2782E+00
10 31 1.4920E-23 3.7171E+00 30 31 5.2667E-12 4.1702E-01
10 32 2.8466E-24 3.9276E+00 30 32 2.5260E-12 5.0965E-01
10 33 4.3736E-25 4.1633E+00 30 33 5.3392E-13 7.1664E-01
10 34 5.7349E-26 4.4148E+00 30 34 6.6950E-14 9.8958E-01
10 35 4.2982E-27 4.7500E+00 30 35 2.9859E-15 1.4140E+00
10 36 7.3294E-28 4.9520E+00 30 36 1.8326E-16 1.7844E+00
10 37 3.7390E-29 5.3420E+00 30 37 1.1507E-17 2.1427E+00
10 38 2.5071E-30 5.6856E+00 30 38 1.5560E-19 2.7166E+00
10 39 2.4084E-31 5.9754E+00 30 39 1.3872E-20 2.9994E+00
10 40 3.5971E-32 6.1965E+00 30 40 5.8839E-22 3.3853E+00
10 41 3.9927E-33 6.4558E+00 30 41 1.1473E-23 3.8789E+00
11 12 1.1117E-14 1.0162E+00 31 32 8.7049E-12 3.3988E-01
11 13 8.0088E-15 1.0725E+00 31 33 3.0925E-12 4.7493E-01
11 14 5.0883E-15 1.1366E+00 31 34 4.2354E-13 7.4367E-01
11 15 2.7065E-15 1.2303E+00 31 35 3.9852E-14 1.0595E+00
11 16 1.8217E-15 1.2821E+00 31 36 1.5830E-15 1.4930E+00
11 17 5.9268E-16 1.4445E+00 31 37 7.1995E-17 1.9075E+00
11 18 2.3725E-16 1.5715E+00 31 38 9.1387E-19 2.4909E+00
11 19 1.0722E-16 1.6799E+00 31 39 4.6151E-20 2.8671E+00
11 20 3.5228E-17 1.8367E+00 31 40 1.7841E-21 3.2648E+00
11 21 1.6589E-17 1.9369E+00 31 41 1.0211E-22 3.5899E+00
11 22 5.2675E-18 2.0942E+00 32 33 6.0387E-12 3.9129E-01
11 23 1.5489E-18 2.2551E+00 32 34 1.3496E-12 5.9144E-01
11 24 6.3090E-19 2.3736E+00 32 35 1.2032E-13 9.1268E-01
11 25 1.8842E-19 2.5308E+00 32 36 6.0680E-15 1.3091E+00
11 26 5.6986E-20 2.6857E+00 32 37 2.7880E-16 1.7238E+00
11 27 1.7468E-20 2.8382E+00 32 38 3.4866E-18 2.3106E+00
11 28 4.0154E-21 3.0261E+00 32 39 2.6549E-19 2.6231E+00
11 29 9.5387E-22 3.2014E+00 32 40 7.9110E-21 3.0681E+00
11 30 1.7494E-22 3.4263E+00 32 41 3.1166E-22 3.4571E+00
11 31 3.6391E-23 3.6258E+00 33 34 4.3448E-12 4.3546E-01
11 32 8.2135E-24 3.8104E+00 33 35 6.0457E-13 6.9780E-01
11 33 1.2530E-24 4.0467E+00 33 36 8.0597E-14 9.6076E-01
11 34 8.8089E-26 4.3902E+00 33 37 5.9674E-15 1.2983E+00
11 35 1.2104E-26 4.6258E+00 33 38 8.0219E-17 1.8821E+00
11 36 1.0687E-27 4.9354E+00 33 39 5.1235E-18 2.2280E+00
11 37 8.3387E-29 5.2616E+00 33 40 2.3432E-19 2.6051E+00
11 38 3.7338E-30 5.6634E+00 33 41 7.7054E-21 3.0228E+00
11 39 4.3531E-31 5.9236E+00 34 35 2.8312E-11 1.6671E-01
11 40 5.2683E-32 6.1757E+00 34 36 1.0598E-11 2.8305E-01
11 41 6.0516E-33 6.4254E+00 34 37 1.4256E-12 5.4302E-01
12 13 1.2597E-14 1.0237E+00 34 38 5.0274E-14 9.8294E-01
12 14 7.2626E-15 1.1067E+00 34 39 5.6757E-15 1.2526E+00
12 15 3.8321E-15 1.2024E+00 34 40 3.6314E-16 1.5926E+00
12 16 2.2059E-15 1.2764E+00 34 41 4.1163E-17 1.8278E+00
12 17 9.8355E-16 1.3941E+00 35 36 5.9527E-11 5.7573E-02
12 18 4.2950E-16 1.5088E+00 35 37 1.0050E-11 2.8291E-01
12 19 2.4239E-16 1.5833E+00 35 38 6.9224E-13 6.3428E-01
12 20 9.3377E-17 1.7159E+00 35 39 6.4186E-14 9.3229E-01
12 21 2.9177E-17 1.8816E+00 35 40 1.0960E-14 1.1286E+00
12 22 1.1829E-17 2.0027E+00 35 41 5.4569E-16 1.5004E+00
12 23 3.4735E-18 2.1641E+00 36 37 2.0427E-11 2.0590E-01
12 24 1.4119E-18 2.2832E+00 36 38 1.6455E-12 5.3500E-01
12 25 4.1961E-19 2.4414E+00 36 39 1.3085E-13 8.5181E-01
12 26 1.2612E-19 2.5974E+00 36 40 1.7326E-14 1.0890E+00
12 27 3.8388E-20 2.7511E+00 36 41 9.8513E-16 1.4505E+00
12 28 1.1688E-20 2.8970E+00 37 38 2.1757E-11 1.7739E-01
12 29 3.3793E-21 3.0429E+00 37 39 1.6662E-12 5.1283E-01
12 30 6.3709E-22 3.2576E+00 37 40 9.5575E-14 8.7830E-01
12 31 8.4167E-23 3.5287E+00 37 41 3.8611E-15 1.2834E+00
12 32 1.7076E-23 3.7322E+00 38 39 1.5671E-11 2.0323E-01
12 33 2.8148E-24 3.9554E+00 38 40 8.5687E-13 5.8192E-01
12 34 2.3801E-25 4.2726E+00 38 41 1.5690E-14 1.1116E+00
12 35 3.7481E-26 4.4909E+00 39 40 5.4921E-12 3.4468E-01
12 36 2.0870E-27 4.8666E+00 39 41 3.7968E-13 6.7653E-01
12 37 1.6058E-28 5.1947E+00 40 41 2.5764E-12 4.3609E-01];
ABdata=[ABdata1(:,1:4); ABdata1(:,5:8)];

out=zeros(42,1);
for i=1:length(ABdata(:,1))
    krate=1.5*ABdata(i, 3).*T.^ABdata(i, 4)/1e6;
    kbf=krate.*exp((O2.e_i(1,ABdata(i, 1))+O2.e_i(1,ABdata(i, 2)+1)...
        -O2.e_i(1,ABdata(i, 1)+1)-O2.e_i(1,ABdata(i, 2)))/k/T);
    R=krate*n(ABdata(i, 1)+1)*n(ABdata(i, 2))...
                                -kbf*n(ABdata(i, 1))*n(ABdata(i, 2)+1);
    out(ABdata(i, 1)+1)=out(ABdata(i, 1)+1)-R;
	out(ABdata(i, 2)) = out(ABdata(i, 2)) - R;
    out(ABdata(i, 1)) = out(ABdata(i, 1)) + R;
	out(ABdata(i, 2)+1)=out(ABdata(i, 2)+1)+R;
end
end