function out=Qel(ne, nm, T)
eedf=load('LoKI-master\Code\Output\RNF\eedf.txt');
me=9.1093837015e-31;
mo2=5.313760980388901e-26;
cross=[  0.000000e+0	3.500000e-21
         1.000000e-3	5.500000e-21
         2.000000e-3	4.100000e-21
         3.000000e-3	4.000000e-21
         5.000000e-3	5.000000e-21
         7.000000e-3	5.800000e-21
         8.500000e-3	6.400000e-21
         1.000000e-2	7.000000e-21
         1.500000e-2	8.700000e-21
         2.000000e-2	9.900000e-21
         3.000000e-2	1.240000e-20
         4.000000e-2	1.440000e-20
         5.000000e-2	1.600000e-20
         7.000000e-2	2.100000e-20
         1.000000e-1	2.500000e-20
         2.000000e-1	3.600000e-20
         3.000000e-1	4.500000e-20
         4.000000e-1	5.200000e-20
         5.000000e-1	5.700000e-20
         6.000000e-1	5.900000e-20
         7.000000e-1	6.100000e-20
         8.000000e-1	6.450000e-20
         9.000000e-1	6.750000e-20
         1.000000e+0	7.200000e-20
         1.100000e+0	7.550000e-20
         1.200000e+0	7.900000e-20
         1.300000e+0	7.900000e-20
         1.400000e+0	7.750000e-20
         1.500000e+0	7.600000e-20
         1.600000e+0	7.450000e-20
         1.700000e+0	7.300000e-20
         1.800000e+0	7.100000e-20
         1.900000e+0	6.900000e-20
         2.000000e+0	6.750000e-20
         2.100000e+0	6.600000e-20
         2.200000e+0	6.500000e-20
         2.300000e+0	6.400000e-20
         2.400000e+0	6.300000e-20
         2.500000e+0	6.100000e-20
         2.600000e+0	6.000000e-20
         2.700000e+0	5.900000e-20
         2.800000e+0	5.800000e-20
         2.900000e+0	5.750000e-20
         3.000000e+0	5.700000e-20
         3.100000e+0	5.650000e-20
         3.200000e+0	5.600000e-20
         3.300000e+0	5.500000e-20
         3.400000e+0	5.480000e-20
         3.500000e+0	5.470000e-20
         3.600000e+0	5.460000e-20
         3.700000e+0	5.450000e-20
         3.800000e+0	5.470000e-20
         3.900000e+0	5.480000e-20
         4.000000e+0	5.500000e-20
         4.100000e+0	5.500000e-20
         4.200000e+0	5.550000e-20
         4.300000e+0	5.550000e-20
         4.400000e+0	5.550000e-20
         4.500000e+0	5.550000e-20
         4.600000e+0	5.560000e-20
         4.700000e+0	5.570000e-20
         4.800000e+0	5.580000e-20
         4.900000e+0	5.590000e-20
         5.000000e+0	5.600000e-20
         5.100000e+0	5.650000e-20
         5.200000e+0	5.670000e-20
         5.300000e+0	5.700000e-20
         5.400000e+0	5.750000e-20
         5.500000e+0	5.800000e-20
         5.600000e+0	5.820000e-20
         5.700000e+0	5.850000e-20
         5.800000e+0	5.870000e-20
         5.900000e+0	5.900000e-20
         6.000000e+0	6.000000e-20
         6.100000e+0	6.100000e-20
         6.200000e+0	6.150000e-20
         6.300000e+0	6.200000e-20
         6.400000e+0	6.250000e-20
         6.500000e+0	6.300000e-20
         6.600000e+0	6.350000e-20
         6.700000e+0	6.400000e-20
         6.800000e+0	6.450000e-20
         6.900000e+0	6.500000e-20
         7.000000e+0	6.600000e-20
         7.100000e+0	6.650000e-20
         7.200000e+0	6.700000e-20
         7.300000e+0	6.750000e-20
         7.400000e+0	6.800000e-20
         7.500000e+0	6.850000e-20
         7.600000e+0	6.900000e-20
         7.700000e+0	6.950000e-20
         7.800000e+0	6.970000e-20
         7.900000e+0	7.000000e-20
         8.000000e+0	7.100000e-20
         8.100000e+0	7.150000e-20
         8.200000e+0	7.170000e-20
         8.300000e+0	7.200000e-20
         8.400000e+0	7.300000e-20
         8.500000e+0	7.400000e-20
         8.600000e+0	7.500000e-20
         8.700000e+0	7.550000e-20
         8.800000e+0	7.600000e-20
         8.900000e+0	7.650000e-20
         9.000000e+0	7.700000e-20
         9.100000e+0	7.720000e-20
         9.200000e+0	7.750000e-20
         9.300000e+0	7.800000e-20
         9.400000e+0	7.850000e-20
         9.500000e+0	7.900000e-20
         9.600000e+0	7.950000e-20
         9.700000e+0	7.970000e-20
         9.800000e+0	7.990000e-20
         9.900000e+0	8.000000e-20
         1.000000e+1	8.000000e-20
         1.010000e+1	8.050000e-20
         1.020000e+1	8.100000e-20
         1.030000e+1	8.150000e-20
         1.040000e+1	8.200000e-20
         1.050000e+1	8.250000e-20
         1.060000e+1	8.260000e-20
         1.070000e+1	8.270000e-20
         1.080000e+1	8.280000e-20
         1.090000e+1	8.290000e-20
         1.100000e+1	8.300000e-20
         1.110000e+1	8.330000e-20
         1.120000e+1	8.360000e-20
         1.130000e+1	8.390000e-20
         1.140000e+1	8.420000e-20
         1.150000e+1	8.450000e-20
         1.160000e+1	8.460000e-20
         1.170000e+1	8.470000e-20
         1.180000e+1	8.480000e-20
         1.190000e+1	8.490000e-20
         1.200000e+1	8.500000e-20
         1.210000e+1	8.540000e-20
         1.220000e+1	8.580000e-20
         1.230000e+1	8.620000e-20
         1.240000e+1	8.660000e-20
         1.250000e+1	8.700000e-20
         1.260000e+1	8.710000e-20
         1.270000e+1	8.720000e-20
         1.280000e+1	8.730000e-20
         1.290000e+1	8.740000e-20
         1.300000e+1	8.750000e-20
         1.310000e+1	8.750000e-20
         1.320000e+1	8.750000e-20
         1.330000e+1	8.760000e-20
         1.340000e+1	8.760000e-20
         1.350000e+1	8.760000e-20
         1.360000e+1	8.760000e-20
         1.370000e+1	8.760000e-20
         1.380000e+1	8.770000e-20
         1.390000e+1	8.770000e-20
         1.400000e+1	8.770000e-20
         1.410000e+1	8.770000e-20
         1.420000e+1	8.770000e-20
         1.430000e+1	8.780000e-20
         1.440000e+1	8.780000e-20
         1.450000e+1	8.780000e-20
         1.460000e+1	8.780000e-20
         1.470000e+1	8.790000e-20
         1.480000e+1	8.790000e-20
         1.490000e+1	8.800000e-20
         1.500000e+1	8.800000e-20
         1.510000e+1	8.800000e-20
         1.520000e+1	8.800000e-20
         1.530000e+1	8.800000e-20
         1.540000e+1	8.800000e-20
         1.550000e+1	8.800000e-20
         1.560000e+1	8.800000e-20
         1.570000e+1	8.790000e-20
         1.580000e+1	8.790000e-20
         1.590000e+1	8.780000e-20
         1.600000e+1	8.780000e-20
         1.610000e+1	8.770000e-20
         1.620000e+1	8.770000e-20
         1.630000e+1	8.760000e-20
         1.640000e+1	8.760000e-20
         1.650000e+1	8.750000e-20
         1.660000e+1	8.740000e-20
         1.670000e+1	8.730000e-20
         1.680000e+1	8.720000e-20
         1.690000e+1	8.710000e-20
         1.700000e+1	8.700000e-20
         1.710000e+1	8.690000e-20
         1.720000e+1	8.680000e-20
         1.730000e+1	8.670000e-20
         1.740000e+1	8.660000e-20
         1.750000e+1	8.650000e-20
         1.760000e+1	8.640000e-20
         1.770000e+1	8.630000e-20
         1.780000e+1	8.620000e-20
         1.790000e+1	8.610000e-20
         1.800000e+1	8.600000e-20
         1.810000e+1	8.600000e-20
         1.820000e+1	8.600000e-20
         1.830000e+1	8.600000e-20
         1.840000e+1	8.600000e-20
         1.850000e+1	8.600000e-20
         1.860000e+1	8.590000e-20
         1.870000e+1	8.580000e-20
         1.880000e+1	8.570000e-20
         1.890000e+1	8.560000e-20
         1.900000e+1	8.550000e-20
         1.910000e+1	8.560000e-20
         1.920000e+1	8.570000e-20
         1.930000e+1	8.580000e-20
         1.940000e+1	8.590000e-20
         1.950000e+1	8.600000e-20
         1.960000e+1	8.600000e-20
         1.970000e+1	8.600000e-20
         1.980000e+1	8.600000e-20
         1.990000e+1	8.600000e-20
         2.000000e+1	8.600000e-20
         2.010000e+1	8.580000e-20
         2.020000e+1	8.560000e-20
         2.030000e+1	8.540000e-20
         2.040000e+1	8.520000e-20
         2.050000e+1	8.500000e-20
         2.060000e+1	8.490000e-20
         2.070000e+1	8.480000e-20
         2.080000e+1	8.470000e-20
         2.090000e+1	8.460000e-20
         2.100000e+1	8.450000e-20
         2.110000e+1	8.440000e-20
         2.120000e+1	8.430000e-20
         2.130000e+1	8.420000e-20
         2.140000e+1	8.410000e-20
         2.150000e+1	8.400000e-20
         2.160000e+1	8.390000e-20
         2.170000e+1	8.380000e-20
         2.180000e+1	8.370000e-20
         2.190000e+1	8.360000e-20
         2.200000e+1	8.350000e-20
         2.210000e+1	8.340000e-20
         2.220000e+1	8.330000e-20
         2.230000e+1	8.320000e-20
         2.240000e+1	8.310000e-20
         2.250000e+1	8.300000e-20
         2.260000e+1	8.300000e-20
         2.270000e+1	8.290000e-20
         2.280000e+1	8.290000e-20
         2.290000e+1	8.280000e-20
         2.300000e+1	8.280000e-20
         2.310000e+1	8.280000e-20
         2.320000e+1	8.280000e-20
         2.330000e+1	8.270000e-20
         2.340000e+1	8.270000e-20
         2.350000e+1	8.270000e-20
         2.360000e+1	8.260000e-20
         2.370000e+1	8.250000e-20
         2.380000e+1	8.250000e-20
         2.390000e+1	8.240000e-20
         2.400000e+1	8.230000e-20
         2.410000e+1	8.230000e-20
         2.420000e+1	8.220000e-20
         2.430000e+1	8.220000e-20
         2.440000e+1	8.210000e-20
         2.450000e+1	8.210000e-20
         2.460000e+1	8.210000e-20
         2.470000e+1	8.210000e-20
         2.480000e+1	8.200000e-20
         2.490000e+1	8.200000e-20
         2.500000e+1	8.200000e-20
         2.510000e+1	8.190000e-20
         2.520000e+1	8.180000e-20
         2.530000e+1	8.170000e-20
         2.540000e+1	8.160000e-20
         2.550000e+1	8.150000e-20
         2.560000e+1	8.150000e-20
         2.570000e+1	8.140000e-20
         2.580000e+1	8.140000e-20
         2.590000e+1	8.130000e-20
         2.600000e+1	8.130000e-20
         2.610000e+1	8.130000e-20
         2.620000e+1	8.120000e-20
         2.630000e+1	8.120000e-20
         2.640000e+1	8.110000e-20
         2.650000e+1	8.110000e-20
         2.660000e+1	8.110000e-20
         2.670000e+1	8.100000e-20
         2.680000e+1	8.100000e-20
         2.690000e+1	8.090000e-20
         2.700000e+1	8.090000e-20
         2.710000e+1	8.090000e-20
         2.720000e+1	8.080000e-20
         2.730000e+1	8.080000e-20
         2.740000e+1	8.070000e-20
         2.750000e+1	8.070000e-20
         2.760000e+1	8.070000e-20
         2.770000e+1	8.070000e-20
         2.780000e+1	8.060000e-20
         2.790000e+1	8.060000e-20
         2.800000e+1	8.060000e-20
         2.810000e+1	8.060000e-20
         2.820000e+1	8.060000e-20
         2.830000e+1	8.050000e-20
         2.840000e+1	8.050000e-20
         2.850000e+1	8.050000e-20
         2.860000e+1	8.050000e-20
         2.870000e+1	8.040000e-20
         2.880000e+1	8.040000e-20
         2.890000e+1	8.030000e-20
         2.900000e+1	8.030000e-20
         2.910000e+1	8.030000e-20
         2.920000e+1	8.020000e-20
         2.930000e+1	8.020000e-20
         2.940000e+1	8.010000e-20
         2.950000e+1	8.010000e-20
         2.960000e+1	8.010000e-20
         2.970000e+1	8.010000e-20
         2.980000e+1	8.000000e-20
         2.990000e+1	8.000000e-20
         3.000000e+1	8.000000e-20
         3.010000e+1	7.990000e-20
         3.020000e+1	7.980000e-20
         3.030000e+1	7.970000e-20
         3.040000e+1	7.960000e-20
         3.050000e+1	7.950000e-20
         3.060000e+1	7.950000e-20
         3.070000e+1	7.940000e-20
         3.080000e+1	7.940000e-20
         3.090000e+1	7.930000e-20
         3.100000e+1	7.930000e-20
         3.110000e+1	7.930000e-20
         3.120000e+1	7.920000e-20
         3.130000e+1	7.920000e-20
         3.140000e+1	7.910000e-20
         3.150000e+1	7.910000e-20
         3.160000e+1	7.910000e-20
         3.170000e+1	7.900000e-20
         3.180000e+1	7.900000e-20
         3.190000e+1	7.890000e-20
         3.200000e+1	7.890000e-20
         3.210000e+1	7.890000e-20
         3.220000e+1	7.880000e-20
         3.230000e+1	7.880000e-20
         3.240000e+1	7.870000e-20
         3.250000e+1	7.870000e-20
         3.260000e+1	7.870000e-20
         3.270000e+1	7.860000e-20
         3.280000e+1	7.860000e-20
         3.290000e+1	7.850000e-20
         3.300000e+1	7.850000e-20
         3.310000e+1	7.850000e-20
         3.320000e+1	7.840000e-20
         3.330000e+1	7.840000e-20
         3.340000e+1	7.830000e-20
         3.350000e+1	7.830000e-20
         3.360000e+1	7.830000e-20
         3.370000e+1	7.820000e-20
         3.380000e+1	7.820000e-20
         3.390000e+1	7.810000e-20
         3.400000e+1	7.810000e-20
         3.410000e+1	7.810000e-20
         3.420000e+1	7.800000e-20
         3.430000e+1	7.800000e-20
         3.440000e+1	7.790000e-20
         3.450000e+1	7.790000e-20
         3.460000e+1	7.790000e-20
         3.470000e+1	7.780000e-20
         3.480000e+1	7.780000e-20
         3.490000e+1	7.770000e-20
         3.500000e+1	7.770000e-20
         3.510000e+1	7.770000e-20
         3.520000e+1	7.760000e-20
         3.530000e+1	7.760000e-20
         3.540000e+1	7.750000e-20
         3.550000e+1	7.750000e-20
         3.560000e+1	7.750000e-20
         3.570000e+1	7.740000e-20
         3.580000e+1	7.740000e-20
         3.590000e+1	7.730000e-20
         3.600000e+1	7.730000e-20
         3.610000e+1	7.730000e-20
         3.620000e+1	7.720000e-20
         3.630000e+1	7.720000e-20
         3.640000e+1	7.710000e-20
         3.650000e+1	7.710000e-20
         3.660000e+1	7.710000e-20
         3.670000e+1	7.700000e-20
         3.680000e+1	7.700000e-20
         3.690000e+1	7.690000e-20
         3.700000e+1	7.690000e-20
         3.710000e+1	7.690000e-20
         3.720000e+1	7.690000e-20
         3.730000e+1	7.680000e-20
         3.740000e+1	7.680000e-20
         3.750000e+1	7.680000e-20
         3.760000e+1	7.680000e-20
         3.770000e+1	7.680000e-20
         3.780000e+1	7.670000e-20
         3.790000e+1	7.670000e-20
         3.800000e+1	7.670000e-20
         3.810000e+1	7.670000e-20
         3.820000e+1	7.670000e-20
         3.830000e+1	7.660000e-20
         3.840000e+1	7.660000e-20
         3.850000e+1	7.660000e-20
         3.860000e+1	7.660000e-20
         3.870000e+1	7.660000e-20
         3.880000e+1	7.650000e-20
         3.890000e+1	7.650000e-20
         3.900000e+1	7.650000e-20
         3.910000e+1	7.650000e-20
         3.920000e+1	7.650000e-20
         3.930000e+1	7.650000e-20
         3.940000e+1	7.650000e-20
         3.950000e+1	7.650000e-20
         3.960000e+1	7.650000e-20
         3.970000e+1	7.650000e-20
         3.980000e+1	7.650000e-20
         3.990000e+1	7.650000e-20
         4.000000e+1	7.650000e-20
         5.000000e+1	7.700000e-20
         7.500000e+1	6.800000e-20
         1.000000e+2	6.500000e-20
         1.500000e+2	6.700000e-20
         2.000000e+2	6.000000e-20
         3.000000e+2	4.900000e-20
         5.000000e+2	3.600000e-20
         7.000000e+2	2.900000e-20
         1.000000e+3	2.120000e-20
         1.500000e+3	1.480000e-20
         2.000000e+3	1.140000e-20
         3.000000e+3	7.900000e-21
         5.000000e+3	5.100000e-21
         7.000000e+3	3.800000e-21
         1.000000e+4	2.800000e-21];
cross_int=interp1(cross(:,1), cross(:,2), eedf(:,1));
df=gradient(eedf(:,2), eedf(:,1));
K=1.3807e-23*6.242e+18;
val=eedf(:,1).^2.*cross_int.*(eedf(:,1)+K*T*df);
Q_int=trapz(eedf(:,1), val);
out=ne*nm*sqrt(2/me)/mo2*Q_int;
end