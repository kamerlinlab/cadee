#cat batch*.log | sed "s=/home/beat/Downloads/cadee/=\$CADEE_DIR/=g" | sed "s=/home/beat/global/=\$HOME/global/=g" | sed "s=beat-ThinkPad-X1-Carbon-3rd=localhost=g" | sed "s= - 2017-09-24 = - 170924 =g" | sed "s/^000//g" | sed "s=/home/beat=\$HOME=g" > finalized.log 
