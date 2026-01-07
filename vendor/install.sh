## SNP2HLA

HOME=SNP2HLA
mkdir -p $HOME && cd $HOME

snp2hla=https://software.broadinstitute.org/mpg/snp2hla/data/SNP2HLA_package_v1.0.3.tar.gz
beagle=https://faculty.washington.edu/browning/beagle/recent.versions/beagle_3.0.4_05May09.zip
beagle2linkage=https://faculty.washington.edu/browning/beagle_utilities/beagle2linkage.jar

wget $snp2hla
wget $beagle
wget $beagle2linkage

tar xvfz `basename $snp2hla`
unzip `basename $beagle`

HOME=home
mkdir -p $HOME && cd $HOME

ln -s ../beagle.3.0.4/beagle.jar .
ln -s ../beagle.3.0.4/utility/linkage2beagle.jar .
ln -s ../beagle2linkage.jar .
ln -s ../SNP2HLA_package_v1.0.3/SNP2HLA/* .
ln -s ../SNP2HLA_package_v1.0.3/Pan-Asian/* .

## xHLA
singularity pull xHLA.sif docker://humanlongevity/hla

## OptiType
singularity pull OptiType.sif docker://fred2/optitype

