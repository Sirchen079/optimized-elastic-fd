function ewm_ensure_dir(pathName)
%EWM_ENSURE_DIR 若目录不存在则创建。

if ~exist(pathName, 'dir')
    mkdir(pathName);
end
end
