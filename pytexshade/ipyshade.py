from pytexshade import shade
from IPython.display import Image
import tempfile
from Bio import Entrez
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment

class shadedmsa(object):
    """Class for storing a shaded image and returning it interactively in jupyter notebook"""

    def __init__(self, msa,shading_modes=['similar'],features=[],title='',legend=True, logo=False,hideseqs=False,splitN=20,setends=[],ruler=False,show_seq_names=True,show_seq_length=True,funcgroups=None,rotate=False,threshold=[80,50],resperline=0,margins=None, density=150,debug=False,startnumber=1):
        temp = tempfile.NamedTemporaryFile(suffix='.png')
        if(debug):
            print("tempfile created: ",temp.name)
        shade.shade_aln2png(msa, filename=temp.name,shading_modes=shading_modes, features=features,title=title,legend=legend, logo=logo,hideseqs=hideseqs,splitN=splitN,setends=setends,ruler=ruler,show_seq_names=show_seq_names,show_seq_length=show_seq_length,funcgroups=funcgroups,rotate=rotate,threshold=threshold,resperline=resperline,margins=margins, density=density,debug=debug,startnumber=startnumber)
        self.img=open(temp.name, 'rb').read()

    def _repr_png_(self):
        return self.img

    
class shadedmsa4plot(object):
    """Class for storing a shaded image and returning it interactively in jupyter notebook:
       version that by default produces an image suitable for annotating an axis of a plot
       i.e. the upper part has a clean border without any annotations
    """

    def __init__(self, msa,shading_modes=['chemical_functional'],features=[],title='',legend=False, logo=False,hideseqs=False,splitN=20,setends=[],ruler='bottom',show_seq_names=False,show_seq_length=False,funcgroups=None,rotate=False,threshold=[80,50],resperline=0,margins=None, density=150,debug=False,startnumber=1):
        temp = tempfile.NamedTemporaryFile(suffix='.png')
        #we'll need to switch all positions of features to bottom forcibly
        fn=features.copy()
        
        for f in fn:
            if f.get('position',False):
                if f['position'][-7:-1]!='bottom':
                    f['position']='bottom'
            else:
                f['position']='bottom'
        
        if(debug):
            print("tempfile created: ",temp.name)
        shade.shade_aln2png(msa, filename=temp.name,shading_modes=shading_modes, features=fn,title=title,legend=legend, logo=logo,hideseqs=hideseqs,splitN=splitN,setends=setends,ruler=ruler,show_seq_names=show_seq_names,show_seq_length=show_seq_length,funcgroups=funcgroups,rotate=rotate,threshold=threshold,resperline=resperline,margins=margins, density=density,debug=debug,startnumber=startnumber)
        self.img=open(temp.name, 'rb').read()

    def _repr_png_(self):
        return self.img

    
    
class shadepdbquick(object):
    """Class for quickly shading a PDB chain with feature info obtained via genbank
    (note, we do not currenlty know were the secondary structure info in genbank comes from
    it is not from PDB file records, it is similar but not exactly the same as information in ss.txt on the PDB web
    site which is claulated using DSSP)
    """

    def __init__(self, pdb_chain_id='1KX5_A',entrez_email='info@example.com',feature_types=['all'],shading_modes=['charge_functional'],title='',legend=True, logo=False,hideseqs=False,splitN=20,setends=[],ruler=True,show_seq_names=True,show_seq_length=True,funcgroups=None,rotate=False,threshold=[80,50],resperline=0,margins=None, density=150,debug=False,startnumber=1):
        temp = tempfile.NamedTemporaryFile(suffix='.png')
        
        #Let's get squence an features
        Entrez.email = entrez_email  # Always tell NCBI who you are
        handle = Entrez.efetch(db="protein", id=pdb_chain_id, rettype="gb", retmode="text")
        record = SeqIO.read(handle, "genbank")
        msar=MultipleSeqAlignment([record])
        msar[0].id='PDB_'+pdb_chain_id

        fn=shade.seqfeat2shadefeat(msar,debug=debug,feature_types=feature_types)
        if(debug):
            print("tempfile created: ",temp.name)
        shade.shade_aln2png(msar, filename=temp.name,shading_modes=shading_modes, features=fn,title=title,legend=legend, logo=logo,hideseqs=hideseqs,splitN=splitN,setends=setends,ruler=ruler,show_seq_names=show_seq_names,show_seq_length=show_seq_length,funcgroups=funcgroups,rotate=rotate,threshold=threshold,resperline=resperline,margins=margins, density=density,debug=debug,startnumber=startnumber)
        self.img=open(temp.name, 'rb').read()

    def _repr_png_(self):
        return self.img


# def getshadedimg(msa,shading_modes=['similar'],features=[],title='',legend=True, logo=False,hideseqs=False,splitN=20,setends=[],ruler=False,show_seq_names=True,show_seq_length=True,funcgroups=None,rotate=False,threshold=[80,50],resperline=0,margins=None, density=150,debug=False):
#     """
#     """

#     temp = tempfile.NamedTemporaryFile(suffix='.png')
#     if(debug):
#         print("tempfile created: ",temp.name)
#     shade.shade_aln2png(msa,filename=temp.name,shading_modes=shading_modes,features=features,title=title,legend=legend, logo=logo,hideseqs=hideseqs,splitN=splitN,setends=setends,ruler=ruler,show_seq_names=show_seq_names,show_seq_length=show_seq_length,funcgroups=funcgroups,rotate=rotate,threshold=threshold,resperline=resperline,margins=margins, density=density,debug=debug)
#     return open(temp.name, 'rb').read()