fid=fopen('c1_1.txt','r');
c11=fscanf(fid,'%f');
fclose(fid)
fid=fopen('c1_2.txt','r');
c12=fscanf(fid,'%f');
fclose(fid)
fid=fopen('c2_2.txt','r');
c22=fscanf(fid,'%f');
fclose(fid)
fid=fopen('c2_1.txt','r');
c21=fscanf(fid,'%f');
fclose(fid)

%% Fast Deconvolution beta=0.1
fid=fopen('InvFD1_1.txt','r');
FD11=fscanf(fid,'%f');
fclose(fid);
fid=fopen('InvFD1_2.txt','r');
FD12=fscanf(fid,'%f');
fclose(fid);
fid=fopen('InvFD2_2.txt','r');
FD22=fscanf(fid,'%f');
fclose(fid);
fid=fopen('InvFD2_1.txt','r');
FD21=fscanf(fid,'%f');
fclose(fid);

s11=conv(c11,FD11)+conv(c21,FD12);
figure;plot(s11);
s21=conv(c11,FD21)+conv(c21,FD22);
figure;plot(s21);
s22=conv(c22,FD22)+conv(c12,FD21);
figure;plot(s22);
s12=conv(c22,FD12)+conv(c12,FD11);
figure;plot(s12);


%% LSE
fid=fopen('InvLSE1_1.txt','r');
LSE11=fscanf(fid,'%f');
fclose(fid);
fid=fopen('InvLSE1_2.txt','r');
LSE12=fscanf(fid,'%f');
fclose(fid);
fid=fopen('InvLSE2_2.txt','r');
LSE22=fscanf(fid,'%f');
fclose(fid);
fid=fopen('InvLSE2_1.txt','r');
LSE21=fscanf(fid,'%f');
fclose(fid);

s11=conv(c11,LSE11)+conv(c21,LSE12);
figure;plot(s11);
s21=conv(c11,LSE21)+conv(c21,LSE22);
figure;plot(s21);
s22=conv(c22,LSE22)+conv(c12,LSE21);
figure;plot(s22);
s12=conv(c22,LSE12)+conv(c12,LSE11);
figure;plot(s12);

%% Fast Deconvolution beta=0.9
fid=fopen('InvFD091_1.txt','r');
FD11=fscanf(fid,'%f');
fclose(fid);
fid=fopen('InvFD091_2.txt','r');
FD12=fscanf(fid,'%f');
fclose(fid);
fid=fopen('InvFD092_2.txt','r');
FD22=fscanf(fid,'%f');
fclose(fid);
fid=fopen('InvFD092_1.txt','r');
FD21=fscanf(fid,'%f');
fclose(fid);

s11=conv(c11,FD11)+conv(c21,FD12);
figure;plot(s11);
s21=conv(c11,FD21)+conv(c21,FD22);
figure;plot(s21);
s22=conv(c22,FD22)+conv(c12,FD21);
figure;plot(s22);
s12=conv(c22,FD12)+conv(c12,FD11);
figure;plot(s12);

%% Fast Deconvolution beta=0.9
fid=fopen('InvFD0011_1.txt','r');
FD11=fscanf(fid,'%f');
fclose(fid);
fid=fopen('InvFD0011_2.txt','r');
FD12=fscanf(fid,'%f');
fclose(fid);
fid=fopen('InvFD0012_2.txt','r');
FD22=fscanf(fid,'%f');
fclose(fid);
fid=fopen('InvFD0012_1.txt','r');
FD21=fscanf(fid,'%f');
fclose(fid);

s11=conv(c11,FD11)+conv(c21,FD12);
figure;plot(s11);
s21=conv(c11,FD21)+conv(c21,FD22);
figure;plot(s21);
s22=conv(c22,FD22)+conv(c12,FD21);
figure;plot(s22);
s12=conv(c22,FD12)+conv(c12,FD11);
figure;plot(s12);
