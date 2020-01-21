from pytexshade.shade import *
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os
import sys
import tempfile
import math
from Bio.Align import MultipleSeqAlignment
from Bio import AlignIO
from io import StringIO
import pkg_resources
import pickle



def test_shade_aln2png_two_seqs():
	print("Testing a small two sequence alignment")
	human_h2a_z_core=Seq('SRSQRAGLQFPVGRIHRHLKSRTTSHGRVGATAAVYSAAILEYLTAEVLELAGNASKDLKVKRITPRHLQLAIRGDEELDSLI-KATIAGGGVIPHIHKSLIG')
	xenopus_h2a_core=Seq('TRSSRAGLQFPVGRVHRLLRKGNYAE-RVGAGAPVYLAAVLEYLTAEILELAGNAARDNKKTRIIPRHLQLAVRNDEELNKLLGRVTIAQGGVLPNIQSVLLP')
	# human_h2a_z_core=Seq('SRSQRAGLQFPVGRIHRHLKSRTTSHGRVGATAAVYSAAILEYLTAEVLELAGNASKDLKVKRITPRHLQLAIRGDEELDSLIKATIAGGGVIPHIHKSLIG')
	msa=MultipleSeqAlignment([SeqRecord(xenopus_h2a_core,id='H2A',name='H2A'),SeqRecord(human_h2a_z_core,id='H2AZ',name='H2AZ')])
	# features=get_hist_ss_in_aln_for_shade(msa,below=True)
	features=[{'style':'fill:$\\uparrow$','sel':[5,10],'text':'test'}]
	print(features)
	shade_aln2png(msa,filename='default',shading_modes=['charge_functional'], legend=False, features=features,title='',logo=False,hideseqs=False,splitN=20,setends=[],ruler=False,show_seq_names=False,show_seq_length=False)
	size=os.path.getsize('default.png')
	assert size>10000, "output png filesize too small, looks that nothing was produced"

def test_5z3l():
	DATA_PATH = pkg_resources.resource_filename('pytexshade', 'pytexshade/data/')

	for dirname, dirnames, filenames in os.walk('.'):
# print path to all subdirectories first.
		for subdirname in dirnames:
			print(os.path.join(dirname, subdirname))
# print path to all filenames.
		for filename in filenames:
			print(os.path.join(dirname, filename))

	msa_dict = pickle.load( open(os.path.join(DATA_PATH,"5z3l_msa_dict.p"), "rb" ) )
	features_dict = pickle.load( open(os.path.join(DATA_PATH,"5z3l_features_dict.p"), "rb" ) )
	for s in 'ABCDEFGH':
		msa=msa_dict[s]
		features=features_dict[s]
		print("Testing chain %s in 5z3l test"%s)
		shade_aln2png(msa,\
 filename='default',\
 shading_modes=['similar'],# list of shading modes according to TexShade, supported are "similar", "hydropathy_functional", "chemical_functional", "structure_functional", "charge_functional", "diverse"\
 legend=False,features=features,title='',\
 logo=False, #SeqLogo \
 hideseqs=False,#activate \hidseqs command \ 
 splitN=20, #alignment will be split into splitN batches\
 setends=[], # a section of alignment will be displayed between setends[0] and setends[1]\
 ruler=True, # Add a ruler\
 show_seq_names=True,# Show sequence names\
 funcgroups=None, # Hilight functional groups, Tex code should be inserted example funcgroups="\\funcgroup{xxx}{CT}{White}{Green}{upper}{up} \\funcgroup{yyy}{GA}{White}{Blue}{upper}{up}" \
 show_seq_length=False #Show sequence length \
)
		size=os.path.getsize('default.png')
		assert size>10000, "output png filesize too small, looks that nothing was produced"
