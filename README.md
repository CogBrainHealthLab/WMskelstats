# WMskelstats

## Development branch
### Replaced `rpyANTs` with pre-installed ANTs binaries
- reticulate no longer required.
- ANTs can be downloaded and installed from https://github.com/ANTsX/ANTs/releases
- ANTs is already installed on NTU HPC (`module load ANTs/2.5.2`)
- Ran the following code successfully. Output files @ `/scratch/junhong.yu/HBN_output`

```
dat=qsi_extract(inputdir="/scratch/junhong.yu/HBN",
              outputdir="/scratch/junhong.yu/HBN_output",
              metrics=c('fa', 'md'),
              skeleton_fathreshold=0.2,
              dti_tensor='dtiTensor',
              dti_L=1,
              nthread=16,
              keep_maps=FALSE,
              silent=FALSE)
```

### To install on NTU HPC (WildFly terminal)
I had to run this on the Wildfly terminal because for some reason `qintel_viz` is not connected to the internet. The above command was ran on `qintel_viz` , after the installation in the Wildfly terminal.


Load modules and open an R session with:
```
module load r/4.6.0
module load gnu/gcc-12.3
ANTs/2.5.2
R
```
Note, any version of R<4.5.0 will not work. 

Install `WMskelstats` (development branch) within this R session:
```
remotes::install_github("cogbrainhealthlab/WMskelstats@development")
```
### To do (03/10/2026)
- need to compare output data with those generated using `rpyANTs`. I could not get `rpyANTs` to work on NTU HPC. I suppose they should be identical since we are using the pre-installed ANTs binary instead a custom-made ANTs replacement function.
- previously, input files were not correctly detected when the input dir contains ses-* sub-directories. In this version of `qsi_extract()`, input paths are modified to detect files correctly within ses-* sub-directory, but these paths may not work when there are no ses-* subdirectories.
- update docs: instructions to install ANTs binaries.

