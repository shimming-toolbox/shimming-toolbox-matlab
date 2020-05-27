

Dependencies:

dcm2niix version v1.0.20200331
Note: older versions do not properly convert json files for the phase data
(one element of ImageType is missing).
 
1)
```
matlab -nodisplay -nosplash -r 'experiment;exit'
```



Octave Notes
------------

Dependencies

1) DICOM

```
octave --eval 'pkg install -forge dicom'
echo "pkg load dicom" >> ~/.octaverc
```

2) NIfTI

There's an in-progress implementation; it looks very complete but has yet to become an official octave package: https://savannah.gnu.org/patch/?9853
To install it manually:

2a. JSONLab
```
curl -JLO https://github.com/fangq/jsonlab/archive/v1.9.8.tar.gz
tar -zxvf jsonlab-1.9.8.tar.gz -C ~/octave/
echo "addpath('~/octave/jsonlab-1.9.8')" >> ~/.octaverc
octave --eval 'loadjson'
```

2b.

```
curl -JLO https://github.com/fangq/zmat/archive/v0.9.2.tar.gz
tar -zxvf zmat-0.9.2.tar.gz -C ~/octave
cp ~/octave/zmat-0.9.2/octave/linux64/zipmat.mex ~/octave/zmat-0.9.2/ # compiled, platform-specific extension module; matlab doesn't handle this well so it's exxxxtra manual.
echo "addpath('~/octave/zmat-0.9.2')" >> ~/.octaverc
octave --eval 'zmat'
```

EDIT:

..this stopped working suddenly (for ArchLinux), failing ? at least on ArchLinux. To fix it I followed the compilation instructions on https://github.com/fangq/zmat/releases:

```
git clone https://github.com/fangq/zmat.git zmat
cd zmat/src/
make clean oct
mv ../zipmat.mex  ~/octave/zmat-0.9.2/
```


2c
```
curl -JLO https://github.com/fangq/jnifti/archive/jnifti_toobox_v0.5.tar.gz
# this is sort of a weird way to install this, but this thing isn't packaged conventionally
# and it actually is two packages, a "matlab" package (which is matlab/octave) implementing JSON-NIfTI,
# and an "octave" package implemented backwards compatibility with Matlab's NIfTI library package.
mkdir -p ~/octave/jnifti && tar -C ~/octave/jnifti -zxvf jnifti-jnifti_toobox_v0.5.tar.gz jnifti-jnifti_toobox_v0.5/lib/ --strip-components 2
echo "addpath(genpath('~/octave/jnifti'))" >> ~/.octaverc
octave --eval 'niftiinfo'
```


3 Image

```
octave --eval 'pkg install -forge image'
echo "pkg load image" >> ~/.octaverc
octave --eval 'mat2gray'
```
