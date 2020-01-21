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

def test_shade_aln2png():
	human_h2a_z_core=Seq('SRSQRAGLQFPVGRIHRHLKSRTTSHGRVGATAAVYSAAILEYLTAEVLELAGNASKDLKVKRITPRHLQLAIRGDEELDSLI-KATIAGGGVIPHIHKSLIG')
	xenopus_h2a_core=Seq('TRSSRAGLQFPVGRVHRLLRKGNYAE-RVGAGAPVYLAAVLEYLTAEILELAGNAARDNKKTRIIPRHLQLAVRNDEELNKLLGRVTIAQGGVLPNIQSVLLP')
    	# human_h2a_z_core=Seq('SRSQRAGLQFPVGRIHRHLKSRTTSHGRVGATAAVYSAAILEYLTAEVLELAGNASKDLKVKRITPRHLQLAIRGDEELDSLIKATIAGGGVIPHIHKSLIG')
	msa=MultipleSeqAlignment([SeqRecord(xenopus_h2a_core,id='H2A',name='H2A'),SeqRecord(human_h2a_z_core,id='H2AZ',name='H2AZ')])
	# features=get_hist_ss_in_aln_for_shade(msa,below=True)
	features=[{'style':'fill:$\\uparrow$','sel':[5,10],'text':'test'}]
	print(features)
	shade_aln2png(msa,filename='default',shading_modes=['charge_functional'], legend=False, features=features,title='',logo=False,hideseqs=False,splitN=20,setends=[],ruler=False,show_seq_names=False,show_seq_length=False)
    
