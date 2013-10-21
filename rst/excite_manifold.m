clear all;
data=[
-70 11.827
-68.3673 13.4596
-66.7347 15.0922
-65.102 16.7247
-63.4694 18.3573
-61.8367 19.9899
-60.2041 21.6225
-58.5714 23.255
-56.9388 24.8876
-55.3061 26.5201
-53.6735 28.1526
-52.0408 29.7852
-50.4082 31.4177
-48.7755 33.0502
-47.1429 34.6827
-45.5102 36.3152
-43.8776 37.9477
-42.2449 39.5801
-40.6122 41.2125
-38.9796 42.8449
-37.3469 44.4772
-35.7143 46.1095
-34.0816 47.7414
-32.449 49.3714
-30.8163 84.3096
-29.1837 85.1836
-27.551 85.6514
-25.9184 85.9936
-24.2857 86.2927
-22.6531 86.6214
-21.0204 86.9003
-19.3878 87.1993
-17.7551 87.5152
-16.1224 87.8148
-14.4898 88.1522
-12.8571 88.4702
-11.2245 88.8
-9.59184 89.1143
-7.95918 89.5062
-6.32653 89.8427
-4.69388 90.2148
-3.06122 90.6294
-1.42857 90.9739
0.204082 91.3113
1.83673 91.6992
3.46939 92.1536
5.10204 92.5268
6.73469 92.8841
8.36735 93.2618
10 93.6116
];
plot(data(:,1),data(:,2),'bo-');
xlabel('V_{ini}, mV');
ylabel('Response Amplitude, mV');
