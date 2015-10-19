clear all

MFILE_ROOT='/net/fs01/data/wzhang/wenchuan/mfilesfinal';
path([MFILE_ROOT '/fun-spool'],path);
path([MFILE_ROOT '/saclab'],path);
path([MFILE_ROOT '/fileexchange'],path);

% ------- coord conversion and directories parameters --------
%conf_crust2_USGS_gtopo30_t12s34f5_smallfault
%fnm_out='crust.USGS.gtopo30.t12s34f5.synthetic.nearest.nc';

conf_700x300_gtopo30_sed1_crust2_ak135_USGS2x2smooth
fnm_out='seismo.recv.nearest.nc';
fnm_mat='seismo.recv.nearest.mat';
id=1; n1=1; n2=14000; dn=1;

% -------------------- load coord --------------------------
subs0=[1,1,1];subc0=[-1,-1,1];subt0=[1,1,1];
[snapinfo]=locate_snap(fnm_conf,id,'start',subs0,'count',subc0,'stride',subt0);
[x,y,z]=gather_coord(snapinfo,'coorddir',pnm_coord);
nx=size(x,1);ny=size(x,2);nz=size(x,3);

% ----------------- station info -----------------------
recv=[];
recv{end+1}={ 'Chengdu', 30+39/60+54/60/60, 104+4/60+19/60/60 };
recv{end+1}={ 'Mianyang',31+28/60+35/60/60,104+43/60+34/60/60};
recv{end+1}={ 'Deyang', 31+8/60+52/60/60, 104+22/60+31/60/60};
recv{end+1}={ 'Mianzhu',31+20/60+44/60/60,104+11/60+40/60/60};
recv{end+1}={ 'Wenchuan',31+29/60+1/60/60,103+35/60+16/60/60};
recv{end+1}={ 'Beichuan',31.8, 104.45};
recv{end+1}={ 'Qingchuan', 32.58  105.23};
recv{end+1}={ 'Langzhong', 31.58  105.96};

%recv{end+1}={ 'ChongQing',29+30/60+54/60/60,106+32/60+40/60/60};
%recv{end+1}={ 'CD2',   3.091000e+01,  1.037580e+02 }; 
%recv{end+1}={ 'LZH',   3.608700e+01,  1.038440e+02 };
%recv{end+1}={ 'GYA',   2.645900e+01,  1.066640e+02 };
%recv{end+1}={ 'XAN',   3.403900e+01,  1.089210e+02 };
%recv{end+1}={ 'KMI',   2.514800e+01,  1.027470e+02 };
%recv{end+1}={ 'GTA',   3.941100e+01,  9.981400e+01 };
%recv{end+1}={ 'LSA',   2.970000e+01,  9.115000e+01 };
%recv{end+1}={ 'TIY',   3.843000e+01,  1.130170e+02 };
%recv{end+1}={ 'HHC',   4.084900e+01,  1.115640e+02 };
%recv{end+1}={ 'GZH',   2.365000e+01,  1.136500e+02 };
%recv{end+1}={ 'TIA',   3.621100e+01,  1.171240e+02 };
%recv{end+1}={ 'NJ2',   3.205200e+01,  1.188540e+02 };
%recv{end+1}={ 'QZN',   1.902900e+01,  1.098440e+02 };
%recv{end+1}={ 'BJI',   4.001800e+01,  1.161700e+02 };
%recv{end+1}={ 'QZH',   2.494300e+01,  1.185920e+02 };
%recv{end+1}={ 'SSE',   3.109600e+01,  1.211870e+02 };
%recv{end+1}={ 'DL2',   3.890600e+01,  1.216280e+02 };
%recv{end+1}={ 'WMQ',   4.381400e+01,  8.770500e+01 };
%recv{end+1}={ 'SNY',   4.182800e+01,  1.235780e+02 };
%recv{end+1}={ 'CN2',   4.380100e+01,  1.254480e+02 };
%recv{end+1}={ 'MDJ',   4.461600e+01,  1.295920e+02 };
%recv{end+1}={ 'ERIC',   32.4592     ,  96.394322    };

nrecv=length(recv);
for n=1:nrecv
    NKNMLEN=length(recv{n}{1});
    lat(n)=recv{n}{2}; lon(n)=recv{n}{3};
    [STX,STY]=geo2cart(lat(n),lon(n),lat0,lon0,alpha_dg,x0,y0);
    [ir,ic]=find(((x-STX).^2+(y-STY).^2)==min(min((x-STX).^2+(y-STY).^2)));
    i=ir(1); j=ic(1);
    indxs{n}=[i,j,1,n1]; indxc{n}=[1,1,1,round((n2-n1+1)/dn)]; indxt{n}=[1,1,1,dn];
    stanm(1:NKNMLEN,n)=recv{n}{1};
    xrecv(n)=STX; yrecv(n)=STY;
end
NKNMLEN=size(stanm,1);

%disp('press enter to continue')
%pause

% --------------------- load data --------------------------
for n=1:nrecv
    disp([ '  retrieve ' num2str(n) 'th recv of ',stanm(:,n)']);

    [snapinfo]=locate_snap(fnm_conf,id,'start',indxs{n},'count',indxc{n},'stride',indxt{n});
    [Vx  ]=retrieve_seismo_snap(snapinfo,id,'Vx','outdir',pnm_out);
    [Vy  ]=retrieve_seismo_snap(snapinfo,id,'Vy','outdir',pnm_out);
    [Vz,T]=retrieve_seismo_snap(snapinfo,id,'Vz','outdir',pnm_out);

    NT=length(T);
    Sx(1:NT,n)=Vx;
    Sy(1:NT,n)=Vy;
    Sz(1:NT,n)=Vz;
end

save(fnm_mat);

% ------------------- save nc -----------------------------
if 1
%my_mode = bitor ( nc_clobber_mode, nc_64bit_offset_mode );
%nc_create_empty ( fnm_out, my_mode );
nc_create_empty(fnm_out);
nc_add_dimension(fnm_out,'number_of_station',nrecv);
nc_add_dimension(fnm_out,'time',NT);
nc_add_dimension(fnm_out,'maxchar',NKNMLEN);

var.Nctype='float';var.Attribute=[];
var.Name='time';var.Dimension={'time'};nc_addvar(fnm_out,var);

var.Name='x';var.Dimension={'number_of_station'};nc_addvar(fnm_out,var);
var.Name='y';var.Dimension={'number_of_station'};nc_addvar(fnm_out,var);
var.Name='latitude';var.Dimension={'number_of_station'};nc_addvar(fnm_out,var);
var.Name='longitude';var.Dimension={'number_of_station'};nc_addvar(fnm_out,var);

var.Name='Vx';var.Dimension={'number_of_station','time'};nc_addvar(fnm_out,var);
var.Name='Vy';var.Dimension={'number_of_station','time'};nc_addvar(fnm_out,var);
var.Name='Vz';var.Dimension={'number_of_station','time'};nc_addvar(fnm_out,var);

var.Nctype='char';var.Attribute=[]; var.Name='station_name';
var.Dimension={'number_of_station','maxchar'};nc_addvar(fnm_out,var);

nc_varput(fnm_out,'time',T);
nc_varput(fnm_out,'x',xrecv);
nc_varput(fnm_out,'y',yrecv);
nc_varput(fnm_out,'latitude',lat);
nc_varput(fnm_out,'longitude',lon);

%nc_varput(fnm_out,'station_name',permute(stanm,[2 1]));
stanm=strjust(stanm','right');
nc_varput(fnm_out,'station_name',stanm);

nc_varput(fnm_out,'Vx',permute(Sx,[2 1]));
nc_varput(fnm_out,'Vy',permute(Sy,[2 1]));
nc_varput(fnm_out,'Vz',permute(Sz,[2 1]));

disp('finished exporting nearest nc')
end

