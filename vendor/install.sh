## SNP2HLA

vendor_dir=$(pwd)

if [ ! -d SNP2HLA ]; then

mkdir -p SNP2HLA && cd SNP2HLA
snp2hla=https://software.broadinstitute.org/mpg/snp2hla/data/SNP2HLA_package_v1.0.3.tar.gz
beagle=https://faculty.washington.edu/browning/beagle/recent.versions/beagle_3.0.4_05May09.zip
beagle2linkage=https://faculty.washington.edu/browning/beagle_utilities/beagle2linkage.jar

wget $snp2hla
wget $beagle
wget $beagle2linkage

tar xvfz `basename $snp2hla`
unzip `basename $beagle`

mkdir -p home && cd home

ln -s ../beagle.3.0.4/beagle.jar .
ln -s ../beagle.3.0.4/utility/linkage2beagle.jar .
ln -s ../beagle2linkage.jar .
ln -s ../SNP2HLA_package_v1.0.3/SNP2HLA/* .
ln -s ../SNP2HLA_package_v1.0.3/Pan-Asian/* .

fi

## DEEP-HLA

cd $vendor_dir

if [ ! -d DEEP-HLA ]; then

mkdir -p DEEP-HLA
git clone https://github.com/tatsuhikonaito/DEEP-HLA.git
cd DEEP-HLA
rm -rf .git

cat <<EOF > environment.yml
name: DEEP-HLA
channels:
  - conda-forge
  - bioconda
  - defaults
dependencies:
  - python=3.7.4
  - setuptools=59.8.0
  - numpy=1.17.2
  - pandas=0.25.1
  - scipy=1.3.1
  - tqdm=4.67.1
  - pip
  - pip:
      - torch==1.4.0
EOF
conda env create -f environment.yml

fi

## xHLA
#singularity pull xHLA.sif docker://humanlongevity/hla

## OptiType
#singularity pull OptiType.sif docker://fred2/optitype

