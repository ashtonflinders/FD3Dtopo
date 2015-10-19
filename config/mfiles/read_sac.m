function read_sac(in,filein)
%script to read sac format file
%pass in filein (filename) and file number (in)
%
global trace stla stlo stel amark 
global tt xx yy
global user0 user1 user2 user3 user4 user5 user6 user7 user8 user9
global kstnm baz dist
global T0 T1 T2 T3 T4 T5 T6 T7 T8 T9
global evla evlo evel evdp

%open file
   [fid,message] = fopen(filein,'r','native');
%    [fid,message] = fopen(filein,'r','b');

%read 70 32 bit floats
   sizef = 70;
   [ff,count] = fread(fid,sizef,'float32');

   delta(in)  = ff(1);
   depmin(in) = ff(2);
   depmax(in) = ff(3);
   scale(in)  = ff(4);
   odelta(in) = ff(5);
   B(in)      = ff(6);
   E(in)      = ff(7);
   O(in)      = ff(8);
   A(in)      = ff(10);
   amark(in)  = ff(9);
   T0(in)     = ff(11);
   T1(in)     = ff(12);
   T2(in)     = ff(13);
   T3(in)     = ff(14);
   T4(in)     = ff(15);
   T5(in)     = ff(16);
   T6(in)     = ff(17);
   T7(in)     = ff(18);
   T8(in)     = ff(19);
   T9(in)     = ff(20);
   F(in)      = ff(21);
   resp0(in)  = ff(22);
   stla(in)   = ff(32);
   stlo(in)   = ff(33);
   stel(in)   = ff(34);
   evla(in)   = ff(36);
   evlo(in)   = ff(37);
   evel(in)   = ff(38);
   evdp(in)   = ff(39);
   user0(in)  = ff(41);
   user1(in)  = ff(42);
   user2(in)  = ff(43);
   user3(in)  = ff(44);
   user4(in)  = ff(45);
   user5(in)  = ff(46);
   user6(in)  = ff(47);
   user7(in)  = ff(48);
   user8(in)  = ff(49);
   user9(in)  = ff(50);
%   enum(in)   = ff(44);
%   tnum(in)   = ff(45);
%  sens(in)   = ff(46);
%   das(in)    = ff(47);
   dist(in)   = ff(51);
   az(in)     = ff(52);
   baz(in)    = ff(53);
   gcarc(in)  = ff(54);
   cmpaz(in)  = ff(58);
   cmpinc(in) = ff(59);
   
%read 35 16 bit integers and logicals
   sizef = 40;
   [ii,count] = fread(fid,sizef,'int');
   nzyear(in) = ii(1);
   nzjday(in) = ii(2);
   nzhour(in) = ii(3);
   nzmin(in)  = ii(4);
   nzsec(in)  = ii(5);
   nzmsec(in) = ii(6);
   npts(in)   = ii(10);
   
%read char
   sizef = 192;
   [cc,count] = fread(fid,sizef,'uchar');
   
   kstnm(in,1:8)  = setstr( cc(1:8) )';
   kevnm(in,1:16) = setstr( cc(9:24) )';
   k0(in,1:6)     = setstr( cc(49:54) )';
   k1(in,1:8)     = setstr( cc(55:62) )';
   k2(in,1:8)     = setstr( cc(64:71) )';
   k3(in,1:8)     = setstr( cc(72:79) )';
   k4(in,1:8)     = setstr( cc(80:87) )';
   k5(in,1:8)     = setstr( cc(88:95) )';
   k6(in,1:8)     = setstr( cc(96:103) )';
   k7(in,1:8)     = setstr( cc(104:111) )';
   k8(in,1:8)     = setstr( cc(112:119) )';
   k9(in,1:8)     = setstr( cc(120:127) )';
    
%read floating 32 bit binary trace data
   sizef = npts(in);
   [data,count] = fread(fid,sizef,'float');
   [mpt,ntrace]=size(data);
  
   tt = B(in) + [0:delta(in):(mpt-1)*delta(in)];
   trace(1:mpt,in) = data;
%   xx(in)=user6(in);
%   yy(in)=user7(in);
%   noise(in)=user9(in);

   [mt,nt]=size(trace);

   flag = 1;

%   plot(tt,trace(:,in),'r-')
%hold on
   
% 05/09/97, Y.shen
    fclose(fid);
%
end
