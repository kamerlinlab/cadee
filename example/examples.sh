i=1

TESTS=1
EDU=1
SCRIPTS=1
DELETE_FILES=0

set -e

export CADEE_DIR="$HOME/Downloads/cadee"

if [[ $DELETE_FILES -eq 1 ]]
then
    echo "Removing old files..."
    /bin/rm -fr "$CADEE_DIR"
    /bin/rm -fr "$HOME/global/cadee_tutorial"
    /bin/rm -fr "$HOME/global/cadee_tutorial_wallclock"
    /bin/rm -fr "$HOME/global/pedagogical_example"
fi

function increment(){
	echo -e "==========================================================================
\n\n\n\n

TIME: $SECONDS
___________________________________________________________________________
NEXT SNIPPET: $i
___________________________________________________________________________
	\n\n\n\n\n
=========================================================================="
	let i+=1
}

if [ $TESTS -eq 1 ]
then

increment

# Code Input Snippet (1)
/bin/bash --version

increment


# Code Input Snippet (2)
sudo apt-get install gfortran openmpi-bin git openbabel
sudo apt-get install mpich gcc python2.7 python-pip 

increment

# Code Input Snippet (3)
cd $HOME
mkdir -p Downloads
cd $HOME/Downloads
git clone https://github.com/kamerlinlab/cadee cadee
cd $HOME/Downloads/cadee 
export CADEE_DIR="$PWD"


git checkout 0.8.5
increment

# Code Input Snippet (4)
python setup.py install --user



increment

# Code Input Snippet (5)
cadee --help

increment

# Code Input Snippet (6)
mkdir testing_example 
cd testing_example
cp -r $CADEE_DIR/example/* .
cadee prep wt.pdb wt.fep wt.qpinp ./libraries/ --template $CADEE_DIR/simpack_templates/simpack_template_0.05ns_15ps_2.5ps_32.5ps.tar.bz2


increment

# Code Input Snippet (7)
set +e
cadee prep wt.pdb wt.fep wt.new.qpinp ./libraries/
if [[ $? -eq 0 ]]
then
    echo "ERROR with SNIPPET 7: Expected Non-Zero Exit-Code."
    exit
fi
set -e


increment

# Code Input Snippet (8)
mkdir -p $HOME/global/cadee_tutorial
cp -r $CADEE_DIR/testing_example/* $HOME/global/cadee_tutorial
# you may need to adjust mpirun
mpirun.mpich -np 5 cadee dyn $HOME/global/cadee_tutorial/wt | tee cadee.log     

i=8

increment

# Code Input Snippet (9)
mkdir -p $HOME/global/cadee_tutorial_wallclock
cp -r $CADEE_DIR/testing_example/* $HOME/global/cadee_tutorial_wallclock
$CADEE_DIR/cadee/tools/pcadee.sh $HOME/global/cadee_tutorial_wallclock/wt 

fi


i=9
if [ $EDU -eq 1 ]
then

increment


# Code Input Snippet (10)
cp -r $CADEE_DIR/example $CADEE_DIR/pedagogical_example
cd $CADEE_DIR/pedagogical_example
cadee prep wt.pdb wt.fep wt.qpinp libraries --alascan --nummuts 48



increment

# Code Input Snippet (11)
cp -r $CADEE_DIR/pedagogical_example $HOME/global
set +e
mpirun.mpich -np 5  cadee dyn $HOME/global/pedagogical_example/ala_scan
set -e

increment

# Code Input Snippet (12)
mpirun.mpich -n 2 cadee dyn $PWD --hij 60.0 --alpha 229.0 --force


increment

# Code Input Snippet (13)
cadee ana alanize cadee.db
firefox index.html 


increment

# Code Input Snippet (14)
cadee ana csv cadee.db activation_barriers.csv	#dG*
cadee ana csv_exo cadee.db free_energy.csv	#ddG
/bin/ls

increment

# Code Input Snippet (15)
set +e
cadee ana cat cadee.db cadee.db1
set -e
/bin/ls *.db

increment

# Code Input Snippet (16)
cd $CADEE_DIR/pedagogical_example
cadee prep wt.pdb wt.fep wt.new.qpinp libraries --libmut 92:SATURATE --libmut 163:SATURATE --libmut 171:SATURATE

increment

# Code Input Snippet (17)
mv libmut point_saturation
cadee prep wt.pdb wt.fep wt.new.qpinp libraries --libmut 92:'AGH' 163:'SPHECA' 171:'WSRLD'

fi


i=17

if [[ $SCRIPTS -eq 1 ]]
then

set +e

increment

# Code Input Snippet (18)
ok=1
msg="\n\n\n\n\n"
python -c 'import cadee' 
if [[ $? -ne 0 ]]
then
    msg="$msg\nUnable to load cadee module!"
    ok=0
fi
cadee > /dev/null
if [[ $? -eq 127 ]]
then
    msg="$msg\nUnable to locate 'cadee' in \$PATH."
    ok=0
fi
if [[ $ok -eq 1 ]]
then
    msg="$msg\nCADEE is installed."
else
    msg="$msg\nCADEE is *NOT* installed!
DIAGNOSIS:"
fi 
echo -e "\n\n$msg"

increment

# Code Input Snippet (19)
# ONLY FOR DEBIAN/UBUNTU
echo $PATH | grep "$HOME/.local/bin" || echo 'Please fix $PATH.'

increment

# Code Input Snippet (20)
# Ubuntu / Debian ONLY
if [ -d "$HOME/.local/bin" ] ; then
    ok=0
    echo $PATH | grep -q "$HOME/.local/bin:" && ok=1
    echo $PATH | grep -q ":$HOME/.local/bin" && ok=1
    if [ $ok -eq 0 ]
    then
        echo 'export PATH="$HOME/.local/bin:$PATH"' >> $HOME/.profile
        source $HOME/.profile
        echo 'Added $HOME/.local/bin to .profile.'
    else
        echo 'Stop: $HOME/.local/bin is already in your $PATH.'
    fi
else
    echo 'Stop: $HOME/.local/bin is not a directory. '
fi

increment

# Code Input Snippet (20)
# Ubuntu / Debian ONLY
if [ -d "$HOME/.local/bin" ] ; then
    ok=0
    echo $PATH | grep -q "$HOME/.local/bin:" && ok=1
    echo $PATH | grep -q ":$HOME/.local/bin" && ok=1
    if [ $ok -eq 0 ]
    then
        echo 'export PATH="$HOME/.local/bin:$PATH"' >> $HOME/.profile
        source $HOME/.profile
        echo 'Added $HOME/.local/bin to .profile.'
    else
        echo 'Stop: $HOME/.local/bin is already in your $PATH.'
    fi
else
    echo 'Stop: $HOME/.local/bin is not a directory. '
fi

increment

# Code Input Snippet (21)	
cadee tool repair_simpack wt_0.tar

fi
