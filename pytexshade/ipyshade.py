from pytexshade import shade
from IPython.display import Image
import tempfile



def shade(msa,filename='default',shading_modes=['similar'],features=[],title='',legend=True, logo=False,hideseqs=False,splitN=20,setends=[],ruler=False,show_seq_names=True,show_seq_length=True,funcgroups=None,rotate=False,threshold=[80,50],resperline=0,margins=None, density=150):
	"""
	"""

	temp = tempfile.NamedTemporaryFile()

	shade.shade_aln2png(msa,filename=temp.name,shading_modes=shading_modes,features=features,title=title,legend=legend, logo=logo,hideseqs=hideseqs,splitN=splitN,setends=setends,ruler=ruler,show_seq_names=show_seq_names,show_seq_length=show_seq_length,funcgroups=funcgroups,rotate=rotate,threshold=threshold,resperline=resperline,margins=margins, density=density)
	Image(filename=temp.name) 