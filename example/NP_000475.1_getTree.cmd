#NEXUS

set warnReset = no;
set increase = auto;
set datastorage=full;
set criterion=parsimony;
execute /Users/jiayinglai/github_vn/example/NP_000475.1.nex;
hsearch nreps=1000 swap=tbr multrees=no;
filter best;
savetrees file=/Users/jiayinglai/github_vn/example/NP_000475.1_trees.nex format=nexus replace=yes root=yes;
quit warnTsave=no;
