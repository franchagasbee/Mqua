#qualimap v2.3

wget https://bitbucket.org/kokonech/qualimap/downloads/qualimap_v2.3.zip
unzip qualimap_v2.3.zip
cd qualimap_v2.3
ln -s $PWD/qualimap ~/bin

#vt

git clone https://github.com/atks/vt.git
cd vt
git submodule update --init --recursive
make -j 10 
make test
ln -s -f $PWD/vt ~/bin