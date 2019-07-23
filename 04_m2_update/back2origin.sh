grep diff m2.update.diff | awk '{print $3}' | perl -e 'while(<>){chomp; s/.*\/src\/main/src\/main/; /.*\/(.*)/; `cp update_package/original_classes/$1 /home/chenxi/sd/test/4070-b6a630a/$_`;}'
