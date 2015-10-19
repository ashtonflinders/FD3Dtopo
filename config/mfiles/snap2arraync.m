clear all

MFILE_ROOT='/net/fs01/data/wzhang/wenchuan/mfilesfinal';
path([MFILE_ROOT '/fun-spool'],path);
path([MFILE_ROOT '/saclab'],path);
path([MFILE_ROOT '/fileexchange'],path);

% ------- coord conversion and directories parameters --------
conf_700x300_gtopo30_sed1_crust2_ak135_USGS2x2smooth

id=1; 
%subs=[201,1,1,1];subc=[1,-1,1,-1];subt=[1,1,1,1];
%fnm_out='seismo.array.I201.nc';
%subs=[401,1,1,1];subc=[1,-1,1,-1];subt=[1,1,1,1];
%fnm_out='seismo.array.I401.nc';
subs=[601,1,1,1];subc=[1,-1,1,-1];subt=[1,1,1,1];
fnm_out='seismo.array.I601.nc';
%subs=[801,1,1,1];subc=[1,-1,1,-1];subt=[1,1,1,1];
%fnm_out='seismo.array.I801.nc';
%subs=[1001,1,1,1];subc=[1,-1,1,-1];subt=[1,1,1,1];
%fnm_out='seismo.array.I1001.nc';
%subs=[1201,1,1,1];subc=[1,-1,1,-1];subt=[1,1,1,1];
%fnm_out='seismo.array.I1201.nc';

%subs=[1,201,1,1];subc=[1,-1,1,-1];subt=[1,1,1,1];
%fnm_out='seismo.array.J201.nc';
%subs=[1,331,1,1];subc=[1,-1,1,-1];subt=[1,1,1,1];
%fnm_out='seismo.array.J331.nc';

rcd=[];
rcd{end+1}={1,[ 201,1,1,1],[1,-1,1,-1],[1,1,1,1],'seismo.array.I201.nc'};
rcd{end+1}={1,[ 401,1,1,1],[1,-1,1,-1],[1,1,1,1],'seismo.array.I401.nc'};
rcd{end+1}={1,[ 601,1,1,1],[1,-1,1,-1],[1,1,1,1],'seismo.array.I601.nc'};
rcd{end+1}={1,[ 801,1,1,1],[1,-1,1,-1],[1,1,1,1],'seismo.array.I801.nc'};
rcd{end+1}={1,[1001,1,1,1],[1,-1,1,-1],[1,1,1,1],'seismo.array.I1001.nc'};
rcd{end+1}={1,[1201,1,1,1],[1,-1,1,-1],[1,1,1,1],'seismo.array.I1201.nc'};

rcd{end+1}={1,[1,201,1,1],[-1,1,1,-1],[1,1,1,1],'seismo.array.J201.nc'};
rcd{end+1}={1,[1,331,1,1],[-1,1,1,-1],[1,1,1,1],'seismo.array.J331.nc'};

nrcd=numel(rcd);

for n=3:nrcd
   id=rcd{n}{1};
   subs=rcd{n}{2};subc=rcd{n}{3};subt=rcd{n}{4};
   fnm_out=rcd{n}{5}

% -------------------- load coord --------------------------
[snapinfo]=locate_snap(fnm_conf,id,'start',subs,'count',subc,'stride',subt);

[X,Y,Z]=gather_coord(snapinfo,'coorddir',pnm_coord);
x=squeeze(X);y=squeeze(Y);z=squeeze(Z);

[V]=retrieve_seismo_snap(snapinfo,id,'Vx','outdir',pnm_out);
Vx=squeeze(permute(V,[4 1 2 3]));
[V]=retrieve_seismo_snap(snapinfo,id,'Vy','outdir',pnm_out);
Vy=squeeze(permute(V,[4 1 2 3]));
[V,T]=retrieve_seismo_snap(snapinfo,id,'Vz','outdir',pnm_out);
Vz=squeeze(permute(V,[4 1 2 3]));

NT=length(T);
nrecv=size(Vz,2);

% ------------------- save nc -----------------------------
%my_mode = bitor ( nc_clobber_mode, nc_64bit_offset_mode );
%nc_create_empty ( fnm_out, my_mode );
nc_create_empty(fnm_out);
nc_add_dimension(fnm_out,'number_of_station',nrecv);
nc_add_dimension(fnm_out,'time',NT);
%nc_add_dimension(fnm_out,'maxchar',NKNMLEN);

var.Nctype='float';var.Attribute=[];
var.Name='time';var.Dimension={'time'};nc_addvar(fnm_out,var);

var.Name='x';var.Dimension={'number_of_station'};nc_addvar(fnm_out,var);
var.Name='y';var.Dimension={'number_of_station'};nc_addvar(fnm_out,var);
var.Name='z';var.Dimension={'number_of_station'};nc_addvar(fnm_out,var);
%var.Name='latitude';var.Dimension={'number_of_station'};nc_addvar(fnm_out,var);
%var.Name='longitude';var.Dimension={'number_of_station'};nc_addvar(fnm_out,var);

var.Name='Vx';var.Dimension={'number_of_station','time'};nc_addvar(fnm_out,var);
var.Name='Vy';var.Dimension={'number_of_station','time'};nc_addvar(fnm_out,var);
var.Name='Vz';var.Dimension={'number_of_station','time'};nc_addvar(fnm_out,var);

%var.Nctype='char';var.Attribute=[]; var.Name='station_name';
%var.Dimension={'number_of_station','maxchar'};nc_addvar(fnm_out,var);

nc_varput(fnm_out,'time',T);
nc_varput(fnm_out,'x',x);
nc_varput(fnm_out,'y',y);
nc_varput(fnm_out,'z',y);
%nc_varput(fnm_out,'latitude',lat);
%nc_varput(fnm_out,'longitude',lon);
%nc_varput(fnm_out,'station_name',permute(stanm,[2 1]));

nc_varput(fnm_out,'Vx',permute(Vx,[2 1]));
nc_varput(fnm_out,'Vy',permute(Vy,[2 1]));
nc_varput(fnm_out,'Vz',permute(Vz,[2 1]));

disp('finished exporting nearest nc')

end
