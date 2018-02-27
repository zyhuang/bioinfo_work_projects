function empbs ()
{
    scriptname=$1;
    template=~/bin/template.pbs;
    jobname=$(basename ${scriptname%.*});
    logdir=$(pwd -P)/log;
    if [ ! -f $1 ]; then
        cat $template | sed 's|JOB_NAME|'$jobname'|g' > $scriptname;
        chmod +x $1;
        mkdir -p $logdir;
    fi;
    emacs -nw $1
}
