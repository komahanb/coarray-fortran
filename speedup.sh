#!/bin/bash

function run_test(){
    EXE=./a.out

    # for (( i=1440; i>360; i-=18 )); do
    for (( i=8; i>1; i-=2 )); do
        echo -n "[" $i ","
        cafrun -np $i ${EXE} | grep "Model run time" | sed  -e's/Model run time:/ /;s/$/],/;s/seconds//'
    done
    echo -n "[ 1, "
    cafrun -np 1 ${EXE} | grep "Model run time" | sed  -e's/Model run time:/ /;s/$/]/;s/seconds//'
}

cat <<EOF
d=np.array([
EOF
run_test
cat <<EOF
])
EOF
