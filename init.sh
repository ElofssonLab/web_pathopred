#!/bin/bash
# initialize the working folder
exec_cmd(){
    echo "$*"
    eval "$*"
}

rundir=$(python -c "import os; print(os.path.realpath(os.path.dirname('$0')))")

cd $rundir

filelist="
$rundir/db.sqlite3
"
dirlist="
$rundir/proj/pred/static/tmp
$rundir/proj/pred/static/result
$rundir/proj/pred/static/md5
$rundir/proj/pred/static/log
$rundir/proj/pred/static/log/stat
$rundir/proj/pred/static/log/divided
"

echo "setting up file permissions"
platform_info=`python -mplatform |  tr '[:upper:]' '[:lower:]'`
platform=
case $platform_info in 
    *centos*)platform=centos;;
    *redhat*) platform=redhat;;
    *ubuntu*|*debian*)platform=ubuntu;;
    *darwin)platform=osx;;
    *)platform=other;;
esac


case $platform in 
    centos|redhat) user=apache;group=apache;;
    ubuntu) user=www-data;group=www-data;;
    osx) user=_www;group=_www;;
    other)echo Unrecognized plat form $platform_info; exit 1;;
esac

# change folder permission and add user to the apache group
myuser=$(whoami)
sudo usermod -a -G $group $myuser
sudo chgrp $group $rundir
sudo chmod 775 $rundir

for file in $filelist; do
    if [ -f "$file" ];then
        exec_cmd "sudo chown $user:$group $file"
    fi
done

echo "create dirs"
for dir in  $dirlist; do
    if [ ! -d $dir ];then
        exec_cmd "sudo mkdir -p $dir"
    fi
    exec_cmd "sudo chmod 755 $dir"
    exec_cmd "sudo chown -R $user:$group $dir"
done

logfile_submit=$rundir/proj/pred/static/log/submitted_seq.log
if [ ! -f $logfile_submit ];then
    exec_cmd "sudo touch $logfile_submit"
fi
exec_cmd "sudo chmod 644 $logfile_submit"
exec_cmd "sudo chown $user:$group $logfile_submit"

# fix the settings.py
if [ ! -f $rundir/proj/settings.py -a ! -L $rundir/proj/settings.py ];then
    pushd $rundir/proj; ln -s pro_settings.py settings.py; popd;
fi

path_config=$rundir/proj/pred/config
configfilelist="
$path_config/alert_email.txt
$path_config/black_iplist.txt
$path_config/computenode.txt
$path_config/vip_email.txt
$path_config/forward_email.txt
$path_config/auth_iplist.txt
"
for f in $configfilelist; do 
    touch $f
done


# create example result
example_folder_list="
example_oneseq
example_multiseq
"
pushd $rundir/proj/pred/static/result

for item in $example_folder_list; do
    if [ ! -d $item ]; then
        sudo ln -s ../download/examples/$item  .
    fi
done
popd
