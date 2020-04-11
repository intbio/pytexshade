from pytexshade import shade
from IPython.display import Image
import tempfile
from Bio import Entrez
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
import requests
import io

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

    def __init__(self, msa,annotation='bottom',shading_modes=['chemical_functional'],features=[],title='',legend=False, logo=False,hideseqs=False,splitN=20,setends=[],ruler=None,show_seq_names=False,show_seq_length=False,funcgroups=None,rotate=False,threshold=[80,50],resperline=0,margins=None, density=150,debug=False,startnumber=1):
        temp = tempfile.NamedTemporaryFile(suffix='.png')
        #we'll need to switch all positions of features to bottom or top forcibly
        fn=features.copy()
        if ruler is None:
            if annotation=='bottom':
                ruler='bottom'
            if annotation=='top':
                ruler='top'
        if annotation=='top':
            for f in fn:
                if f.get('position',False):
                    if f['position'][-3:]!='top':
                        f['position']='top'
                else:
                    f['position']='bottom'
        else:
            for f in fn:
                if f.get('position',False):
                    if f['position'][-6:]!='bottom':
                        if debug:
                            print('Setting position to bottom, initial position is %s'%f['position'])
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
    """Class for quickly shading a PDB chain with feature info obtained via genbank (this would be the SEQREC record from PDB)
    (note, we do not currenlty know were the secondary structure info in genbank comes from
    it is not from PDB file records, it is similar but not exactly the same as information in ss.txt on the PDB web
    site which is claulated using DSSP)
    """

    def __init__(self, pdb_chain_id='1KX5_A',entrez_email='info@example.com',feature_types=['all'],add_features=[],force_feature_pos=None,shading_modes=['charge_functional'],title='',legend=True, logo=False,hideseqs=False,splitN=20,setends=[],ruler=True,show_seq_names=True,show_seq_length=True,funcgroups=None,rotate=False,threshold=[80,50],resperline=0,margins=None, density=150,debug=False,startnumber=1):
        temp = tempfile.NamedTemporaryFile(suffix='.png')
        
        self.seqrec=False
        #Let's try the NCBI way if it is a protein seq
        try:
            Entrez.email = entrez_email  # Always tell NCBI who you are
            handle = Entrez.efetch(db="protein", id=pdb_chain_id, rettype="gb", retmode="text")
            record = SeqIO.read(handle, "genbank")
            self.seqrec=record
            if(debug):
                print(record)
            msar=MultipleSeqAlignment([self.seqrec])
            msar[0].id='PDB_'+pdb_chain_id
            fn=shade.seqfeat2shadefeat(msar,debug=debug,feature_types=feature_types,force_feature_pos=force_feature_pos)
            if(debug):
                print(fn)
        except Exception as exc:
            pass
#             print(exc)
        
        if(not self.seqrec):
            #Let's get from pdb
            pdbid=pdb_chain_id.split('_')[0]
            chid=pdb_chain_id.split('_')[1]
            seqrec={}
            h=io.StringIO(requests.get('http://www.rcsb.org/pdb/download/viewFastaFiles.do?structureIdList=%s&compressionType=uncompressed'%pdbid).content.decode("utf-8") )
            for record in SeqIO.parse(h,'fasta'):
                seqrec[record.id.split(':')[1].split('|')[0]]=record
#             print(seqrec)
            self.seqrec=seqrec[chid]
            fn=[]
        
        msar=MultipleSeqAlignment([self.seqrec])
        msar[0].id='PDB_'+pdb_chain_id
        fn.extend(add_features)
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