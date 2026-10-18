function path = HydroPol2D_Write_Run_Manifest(case_root,results,output_root)
%HYDROPOL2D_WRITE_RUN_MANIFEST Record the run and immutable inputs.

arguments
    case_root (1,:) char
    results struct
    output_root (1,:) char
end

config_path=fullfile(case_root,'hydropol_config.json');
config=jsondecode(fileread(config_path));

model_root=fileparts(fileparts(mfilename('fullpath')));
[revision_status,revision]=system(sprintf('git -C "%s" rev-parse HEAD',model_root));
if revision_status~=0, revision='unavailable'; end
[dirty_status,dirty]=system(sprintf('git -C "%s" status --porcelain',model_root));
if dirty_status~=0, dirty='unavailable'; end

relative_inputs={ ...
    'hydropol_config.json', ...
    fullfile('processed_inputs','Mesh','hydropol_quadtree_mesh.nc'), ...
    fullfile('processed_inputs','Mesh','hydropol_quadtree_overlap.nc'), ...
    fullfile('processed_inputs','Static','DEM.tif'), ...
    fullfile('processed_inputs','Static','LULC.tif'), ...
    fullfile('processed_inputs','Static','SOIL.tif'), ...
    fullfile('processed_inputs','Static','DTB.tif'), ...
    fullfile('processed_inputs','Static','GW_table.tif'), ...
    fullfile('processed_inputs','Static','RiverbedHydraulicConductivity.tif')};
inputs=struct('path',{},'sha256',{});
for k=1:numel(relative_inputs)
    candidate=fullfile(case_root,relative_inputs{k});
    if ~isfile(candidate), continue; end
    inputs(end+1)=struct('path',relative_inputs{k},'sha256',sha256_file(candidate)); %#ok<AGROW>
end

latest_s=max(double(results.time_s(:)),[],'omitnan');
start_time=datetime(config.simulation.start_utc,'TimeZone','UTC', ...
    'InputFormat','yyyy-MM-dd''T''HH:mm:ssXXX');
latest_time=start_time+seconds(latest_s);
manifest=struct( ...
    'solver','HydroPol2D MATLAB', ...
    'solver_revision',strtrim(revision), ...
    'solver_worktree_dirty',~isempty(strtrim(dirty)), ...
    'configuration',config, ...
    'input_digests',inputs, ...
    'latest_saved_model_time_s',latest_s, ...
    'latest_saved_model_time_utc',char(string(latest_time,'yyyy-MM-dd''T''HH:mm:ssXXX')), ...
    'saved_output_count',numel(results.time_s));

if ~isfolder(output_root), mkdir(output_root); end
path=fullfile(output_root,'run_manifest.json');
fid=fopen(path,'w');
assert(fid>=0,'HydroPol2D:RunManifest','Could not create %s.',path);
cleanup=onCleanup(@() fclose(fid));
fprintf(fid,'%s\n',jsonencode(manifest,PrettyPrint=true));
end

function digest=sha256_file(path)
fid=fopen(path,'r');
assert(fid>=0,'HydroPol2D:RunManifest','Could not read %s.',path);
cleanup=onCleanup(@() fclose(fid));
bytes=fread(fid,inf,'*uint8');
engine=java.security.MessageDigest.getInstance('SHA-256');
engine.update(bytes);
digest=lower(reshape(dec2hex(typecast(engine.digest(),'uint8'),2).',1,[]));
end
