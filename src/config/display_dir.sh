namevar=${1}
pathvar=${2}
printf "%17s" $namevar
printf "%1s" ' '
printf "%s" ${pathvar##*/augeanstables/}
printf "\n"