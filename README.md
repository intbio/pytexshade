# pytexshade
A python wrapper for [TexShade](https://ctan.org/pkg/texshade?lang=en) sequence alignment shader

## Installation via conda
```
conda isntall -c intbio -c conda-forge pytexshade
```


## Usage
### In Jupyter
For Jupyter use see [example.ipynb](example.ipynb) - THIS IS RECOMMENDED
### Stand alone

```
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment

from pytexshade.shade import shade_aln2png

human_h2a_z_core=Seq('SRSQRAGLQFPVGRIHRHLKSRTTSHGRVGATAAVYSAAILEYLTAEVLELAGNASKDLKVKRITPRHLQLAIRGDEELDSLI-KATIAGGGVIPHIHKSLIG')
xenopus_h2a_core=Seq('TRSSRAGLQFPVGRVHRLLRKGNYAE-RVGAGAPVYLAAVLEYLTAEILELAGNAARDNKKTRIIPRHLQLAVRNDEELNKLLGRVTIAQGGVLPNIQSVLLP')

   msa=MultipleSeqAlignment([SeqRecord(xenopus_h2a_core,id='H2A',name='H2A'),SeqRecord(human_h2a_z_core,id='H2AZ',name='H2AZ')])

features=[{'style':'fill:$\\uparrow$','sel':[5,10],'text':'test'}]
shade_aln2png(msa,filename='shaded.png',shading_modes=['charge_functional'], legend=False,features=features,title='',logo=False,hideseqs=False,splitN=20,setends=[],ruler=False,show_seq_names=False,show_seq_length=False)
    
```
Result is
![shaded.png](shaded.png)


Main function is 
```
shade_aln2png(msa,filename='default',shading_modes=['similar'],features=[],title='',legend=True, logo=False,hideseqs=False,splitN=20,setends=[],ruler=False,show_seq_names=True,show_seq_length=True,funcgroups=None,rotate=False,threshold=[80,50],resperline=0,margins=None, density=150)
```

It will convert MultipleSequenceAlignment to a shaded image.
- shading_modes: similar, ... see write_texshade code
- features - a list of dictionaries:
{'style':'shaderegion','shadeblock','frameblock','helix','loop','-->','---','<--',',-,' - all texshade features types of features
+ frameblock + shaderegion + shadeblock
'position':'top','bottom','ttop','bbottom', etc. if no - automatic
'seqref':number of sequence for selection - default consensus
'sel':[begin,end] - region selection, this is in 0 based numbering (as in msa - we override here the TexShade 1-based numbering)
'text':'text'
'color':'Red' this is for frame block
'thickness':1.5 -1.5pt also for frameblock
'rescol':'Black' this is for \shaderegion{ seqref  }{ selection }{ res.col. }{ shad.col. }
'shadcol':'Green' - the same
funcgroup example fg="\\funcgroup{xxx}{CT}{White}{Green}{upper}{up} \\funcgroup{xxx}{GA}{White}{Blue}{upper}{up}"
}


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
