# pytexshade
A python wrapper for [TexShade](https://ctan.org/pkg/texshade?lang=en) sequence alignment shader

## Installation via conda
```
conda create --name pytexshade
source activate pytexshade
conda install -c conda-forge tectonic
conda install -c conda-forge imagemagick
conda isntall -c intbio pytexshade

```




## Usage


### Alternative installation notes
TeX distro for Mac http://www.tug.org/mactex/morepackages.html
and then add TeXShade via
```
sudo tlmgr update --self
sudo tlmgr install texshade
```

- ImageMagic has to be installed - the `convert` command at least should be available systemwide.
The Ghostscript should be installed - the `gs` command at least should be available systemwide.

Via brew:
```
brew install im
brew install gs
```
