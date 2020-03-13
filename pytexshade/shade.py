# -*- coding: utf-8 -*-
"""
This is a library that makes good images of shaded alignments 
through TeXShade.

Input:
1) alignment (might consist of a single sequence).
2) Shading options
3) Features list
Output:
pdf or image (needed for automaitc plotting)


"""

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os
import sys
import tempfile
import math
from Bio.Align import MultipleSeqAlignment
from Bio import AlignIO
from io import StringIO


def seqfeat2shadefeat(msa,seqref=None,idseqref=True, debug=False):
    """
    converts SeqFeature records from genbank for every sequence in msa to our style feature list
    """
    gb_feat_dict_conv={'sheet':'-->'}
    features=[]
    #In texshade [1,1] will be colored as one residue. this is different from biopython where [1,1] selects nothing!
    for m,i in zip(msa,range(len(msa))):
        if(debug):
            print(m)
        for f in m.features:
            if(debug):
                print(f)
            if seqref:
                sr=seqref
            elif idseqref:
                sr=m.id
            else:
                sr=i+1
            if f.type=='SecStr':
                features.append({'style':f.qualifiers['sec_str_type'][0],'seqref':sr,'sel':[f.location.start,f.location.end-1]})   
            if f.type=='motif':
                features.append({'style':'shaderegion','seqref':sr,'sel':[f.location.start,f.location.end-1],'shadcol':f.qualifiers['color']})
            if f.type=='domain':
                features.append({'style':'shaderegion','seqref':sr,'sel':[f.location.start,f.location.end-1],'shadcol':f.qualifiers['color']})
            if f.type=='Region':
                features.append({'style':'---','seqref':sr,'sel':[f.location.start,f.location.end-1],'color':'black','text':f.qualifiers['region_name'][0],'position':'bottom'})
            if f.type=='frameblock':
                features.append({'style':'frameblock','seqref':sr,'sel':[f.location.start,f.location.end-1],'color':f.qualifiers['color']})
            if f.type=='---':
                features.append({'style':'---','seqref':sr,'sel':[f.location.start,f.location.end-1],'color':f.qualifiers['color'],'text':f.qualifiers['text']})

    return features

def shade_aln2png(msa,filename='default',shading_modes=['similar'],features=[],title='',legend=True, logo=False,hideseqs=False,splitN=20,setends=[],ruler=False,show_seq_names=True,show_seq_length=True,funcgroups=None,rotate=False,threshold=[80,50],resperline=0,margins=None, density=150,debug=False,startnumber=1):
    td=tempfile.TemporaryDirectory()
    TEMP_DIR=td.name
    if(debug):
        print("Created temporaty directory: ", TEMP_DIR)
    intf=os.path.join(TEMP_DIR,'tempshade.pdf')
    shade_aln2pdf(msa,intf,shading_modes,features,title,legend, logo,hideseqs,splitN,setends,ruler,show_seq_names,show_seq_length,funcgroups,threshold,resperline=resperline,debug=debug,startnumber=startnumber)
    #let's use imagemagic
    #margins - add on each side margins %
    if margins:
        m=margins
    else:
        m=0
    if(debug):
        print("DEBUG: copying pdf in current folder for analysis")
        cmd="mkdir -p debug && cp %s debug/tempshade.pdf"%intf
        print(cmd)
        os.system(cmd)
        print("Converting PDF to PNG")
    if rotate:
        cmd='convert -density %d '%density+intf+' -trim -bordercolor White -border %.3f%%x0%% -rotate -90 %s'%(m,filename if filename[-3:]=='png' else filename+'.png')
        if(debug):
            print(cmd)
        os.system(cmd)
    else:
        cmd='convert -density %d '%density+intf+' -trim -bordercolor White -border %.3f%%x0%% %s'%(m,filename if filename[-3:]=='png' else filename+'.png')
        if(debug):
            print(cmd)
        os.system(cmd)
    os.remove(intf)


def shade_aln2pdf(msa,filename='default',shading_modes=['similar'],features=[],title='',legend=True, logo=False,hideseqs=False,splitN=20,setends=[],ruler=False,show_seq_names=True,show_seq_length=True,funcgroups=None,threshold=[80,50],resperline=0,debug=False,startnumber=1):
    """
will convert msa to a shaded pdf.
shading_modes: similar, ... see write_texshade code
features - a list of dictionaries:
{'style':'shaderegion','shadeblock','frameblock','helix','loop','-->','---','<--',',-,' - all texshade features types of features
+ frameblock + shaderegion + shadeblock
'position':'top','bottom','ttop','bbottom', etc. if no - automatic
'seqref':number of sequence for selection - default consensus
'setends':  analougs to setends command in TexShade - set the part of sequence to display - numbering as in TexShade, currently this is buggy in TexShade.
'startnumber' - texshade starnumber parameter in \\stratnumber of \\setends - currently is also buggy.
'sel':[begin,end] - region selection, this is in 0 based numbering (as in msa - we override here the TexShade 1-based numbering)
'text':'text'
'color':'Red' this is for frame block
'thickness':1.5 -1.5pt also for frameblock
'rescol':'Black' this is for \\shaderegion{ seqref  }{ selection }{ res.col. }{ shad.col. }
'shadcol':'Green' - the same}
ruler - show numbers every 10 residues. Can be True - show on top, or 'bottom' - show below.
funcgroup example fg="\\funcgroup{xxx}{CT}{White}{Green}{upper}{up} \\funcgroup{xxx}{GA}{White}{Blue}{upper}{up}"
resperline - number of residues per line, if 0 - then equal to the length of the sequence.
startnumber - hike the starting number for ruler by startnumber-1


    """

    td=tempfile.TemporaryDirectory()
    TEMP_DIR=td.name
    if(debug):
        print("Created temporaty directory: ", TEMP_DIR)

    ns='consensus'
    hideseqs_by_name=[]
    if(len(msa)==1):
        msa=msa[:]
        msa.extend([SeqRecord(msa[0].seq,id='dum',name='dum')])
        hideseqs_by_name.append('dum')
    #####if we are splitting the alignment into blocks - get number of blocks

    a_len=len(msa)
    num=int(a_len/splitN)+1
    while ((a_len-(num-1)*splitN)<2):
        splitN=splitN+1
        num=int(a_len/splitN)+1
    if(debug):
        print("Chosen splitting parameters")
        print(a_len, splitN)

    ### rename any seqs with _ to -
    for i in range(len(msa)):
        if(msa[i].id[0].isdigit()):
            raise Exception('Tex cannot handle sequence ids starting with a number, offending id is %s'%msa[i].id)
        msa[i].id=msa[i].id.replace('_','-')
        msa[i].name=msa[i].name.replace('_','-')


    ####iterate over blocks and create alignment fasta files
    for i in range(num):
        t_aln=msa[(i*splitN):((i+1)*splitN)]
        AlignIO.write(t_aln,open(os.path.join(TEMP_DIR,'alignment%d.fasta'%i),'w'), 'fasta')
        if(debug):
            print("DEBUG: copying fasta file in current folder for analysis")
            cmd="mkdir -p debug && cp %s debug/alignment%d.fasta"%(os.path.join(TEMP_DIR,'alignment%d.fasta'%i),i)
            print(cmd)
            os.system(cmd)
    if resperline==0:
        res_per_line=len(msa[0])
    else:
        res_per_line=resperline

    #prepare feature section

    #alias dict
    aliasf={'alpha':'helix','beta':'-->','sheet':'-->','domain':'loop'}


    features_dict={}
    for i in features:
        sr=str(i.get('seqref','consensus')).replace('_','-')
        if(i['style']=='block' or i['style']=='frameblock'):
            features_dict[sr]=features_dict.get(sr,'')+"\\frameblock{%s}{%d..%d}{%s[%.1fpt]}"%(sr,i['sel'][0]+1,i['sel'][1]+1,i.get('color','Red'),i.get('thickness',1.5))
        elif(i['style']=='shaderegion' or i['style']=='shadeblock'):
            features_dict[sr]=features_dict.get(sr,'')+"\\%s{%s}{%d..%d}{%s}{%s}"%(i['style'],sr,i['sel'][0]+1,i['sel'][1]+1,i.get('rescol','Black'),i.get('shadcol','Green'))
        else:
            features_dict[sr]=features_dict.get(sr,'')+"\\feature{%s}{%s}{%d..%d}{%s}{%s}"%(i.get('position','top'),sr,i['sel'][0]+1,i['sel'][1]+1,aliasf.get(i.get('style','loop'),i.get('style','loop')),i.get('text','').replace('_','-'))
    if(debug):
        print(features_dict)
    a=open(os.path.join(TEMP_DIR,'align.tex'),'w')

    a.write(r"""\documentclass[11pt,landscape]{article}
%\documentclass{standalone}
%\usepackage[a0paper]{geometry}
\usepackage{hyperref}
""")

    lmult=math.ceil(float(len(msa[0]))/float(res_per_line))
    ca_len=a_len*lmult
    h=((ca_len/30.*25 + (2.5 if legend else 0.0)) if (ca_len/30.*25 + (2.5 if legend else 0.0) <25.0) else 25)+4
    if(logo):
        h=h+1
    w=(22/200.*res_per_line+2.5)

    if title:
        h+=1
        if w<len(title)*0.4:
            w=len(title)*0.4
    a.write("""
\\usepackage[paperwidth=%fin, paperheight=%fin]{geometry}
        """%(w,h))

    a.write(r"""
\usepackage{texshade}

\begin{document}""")
    a.write("""
\\pagenumbering{gobble}
\\centering
\\Huge{%s}
\\vspace{-0.5in}"""%title)

    for i in range(num-1):
        features_code=features_dict.get('consensus','')
        for s,ns in zip(msa[(i*splitN):((i+1)*splitN)],range((i*splitN),((i+1)*splitN))):
            features_code+=features_dict.get(s.id,'')
            features_code+=features_dict.get(str(ns+1),'')
        write_texshade(a,os.path.join(TEMP_DIR,'alignment%d.fasta'%i) , features_code, res_per_line,False,shading_modes,logo,hideseqs,setends,ruler,numbering_seq='consensus',hide_ns=False,show_seq_names=show_seq_names,show_seq_length=show_seq_length,hideseqs_by_name=hideseqs_by_name,funcgroups=funcgroups,threshold=threshold,startnumber=startnumber)
    
    i=num-1
    features_code=features_dict.get('consensus','')
    for s,ns in zip(msa[(i*splitN):((i+1)*splitN)],range((i*splitN),((i+1)*splitN))):
        features_code+=features_dict.get(s.id,'')
        features_code+=features_dict.get(str(ns+1),'')
    write_texshade(a,os.path.join(TEMP_DIR,'alignment%d.fasta'%(num-1)) , features_code, res_per_line,legend,shading_modes,logo,hideseqs,setends,ruler,numbering_seq='consensus',hide_ns=False,show_seq_names=show_seq_names,show_seq_length=show_seq_length,hideseqs_by_name=hideseqs_by_name,funcgroups=funcgroups,threshold=threshold,startnumber=startnumber)

    a.write(r"""
\end{document} """)
    a.close()

    # command='pdflatex --file-line-error --synctex=1 -output-directory=%s --save-size=10000  %s/align.tex > /dev/null'%(TEMP_DIR,TEMP_DIR)
    command='tectonic %s/align.tex 2>&1'%(TEMP_DIR)

    if(debug):
        print("DEBUG: copying tex file in current folder for analysis")
        cmd="mkdir -p debug && cp %s/align.tex debug/align.tex"%TEMP_DIR
        print(cmd)
        os.system(cmd)
        print('Launcning Tectonic:')
        
        print(command)
    out=os.popen(command).read()
    if(debug):
        print(out)
    command='mv '+os.path.join(TEMP_DIR,'align.pdf')+' %s'%(filename if filename[-3:]=='pdf' else (filename+'.pdf'))
    if(debug):
        print(command)
    os.system(command)
    # for i in range(num):
        # os.remove(os.path.join(TEMP_DIR,'alignment%d.fasta'%i))
    # os.remove(os.path.join(TEMP_DIR,'align.tex'))




def write_texshade(file_handle,aln_fname,features,res_per_line=120,showlegend=True,shading_modes=['similar'],logo=False,hideseqs=False,setends=[],ruler=False,numbering_seq='consensus',hide_ns=False,show_seq_names=True,show_seq_length=True,hideseqs_by_name=[],funcgroups=None,threshold=[80,50],startnumber=1):
    if((not setends) and (startnumber!=1)):
        raise Exception("With startnumber not equal 1, pls, define setends manually!")
    if(setends and (startnumber!=1)):
        startnumber=startnumber+1 #There is some bug in TexShade currently, I think.
        
    for shading in shading_modes:
        shading=str(shading)
        file_handle.write("""
    \\begin{texshade}{%s}
    \\residuesperline*{%d}
  %% \\allowzero
    """%(aln_fname.split('/')[-1],res_per_line)) # we expect that tex file will be in the same dir as alignment files!!! hense split('/')
       


        file_handle.write(r"""
    \seqtype{P}
    \defconsensus{{}}{*}{upper}
    """)
        #a very dirty hack
        if(setends):
            file_handle.write("""
     \\setends[%d]{%s}{%d..%d}
     """%(startnumber,numbering_seq,setends[0],setends[1]))
                    
            
# I have no ideas why whas that needed . @molsim 22 Jan 2020            
#             if(numbering_seq=='consensus'):
#                 file_handle.write("""
#     \\setends[%d]{%s}{%d..%d}
#     """%(setends[0]+startnumber,numbering_seq,setends[0]+setends[0],setends[1]+setends[0]+setends[0]-1))
#             else:
#                     file_handle.write("""
#     \\setends[%d]{%s}{%d..%d}
#     """%(startnumber,numbering_seq,setends[0],setends[1]))
                    
                    
        else:
            file_handle.write("""
    \\startnumber{%s}{%d}
    """%(numbering_seq,startnumber))
           
        if(ruler):
            if(ruler=='bottom'):
                file_handle.write("""
    \\showruler{bottom}{%s}
    """%numbering_seq)
            else:
                file_handle.write("""
    \\showruler{top}{%s}
    """%numbering_seq)
        if(hide_ns):
            file_handle.write("""
    \\hideseq{%s}
    """%numbering_seq)

        file_handle.write(r"""
    \seqtype{P}
    """)
        if(not show_seq_names):
            file_handle.write(r"""
    \hidenames
    """) 
        if(not show_seq_length):
            file_handle.write(r"""
    \hidenumbering
    """) 
        if((shading=='similar')|(shading=='0')):
            file_handle.write(r"""
    \shadingmode{similar}
    \threshold[%d]{%d}
    """%(threshold[0],threshold[1]))
            if(logo):
                file_handle.write(r"""
    \showsequencelogo{top} \showlogoscale{leftright}
    \namesequencelogo{logo}
    """)

        if((shading=='hydropathy_functional')|(shading=='1')):
            file_handle.write(r"""
    \shadingmode[hydropathy]{functional}
    \shadeallresidues
    \threshold[%d]{%d}
    
    """%(threshold[0],threshold[1]))     
            if(logo):
                file_handle.write(r"""
    \showsequencelogo[hydropathy]{top} \showlogoscale{leftright}

    """)
        if((shading=='chemical_functional')|(shading=='2')):
            file_handle.write(r"""

    \shadingmode[chemical]{functional}
    \shadeallresidues
    """)
            if(logo):
                file_handle.write(r"""
    \showsequencelogo[chemical]{top} \showlogoscale{leftright}

    """)
        if((shading=='structure_functional')|(shading=='3')):
            file_handle.write(r"""
    \shadingmode[structure]{functional}
    \shadeallresidues
    """)

        if((shading=='charge_functional')|(shading=='4')):
            file_handle.write(r"""
    \shadingmode[charge]{functional}
    \shadeallresidues
    """)

        if((shading=='diverse')|(shading=='5')):
            file_handle.write(r"""
    \shadingmode{diverse}

    """)
        if(funcgroups):
            file_handle.write(funcgroups)
            
        if(hideseqs):
            file_handle.write(r"""
    \hideseqs

    """)
        for s in hideseqs_by_name:
            file_handle.write("""
    \\hideseq{%s}
    """%s)

        file_handle.write(r"""

    %\setends{consensus}{1..160}
    %\setends{consensus}{1..160}


    %\feature{ttop}{1}{1..160}{bar:conservation}{}
    %\showfeaturestylename{ttop}{conserv.}
    \ttopspace{-\baselineskip}

    %\feature{top}{1}{1..160}{color:charge}{}
    %\showfeaturestylename{top}{charge}

    %\feature{bottom}{1}{1..160}{color:molweight[ColdHot]}{}
    %\showfeaturestylename{bottom}{weight}

    %\bbottomspace{-\baselineskip}
    %\feature{bbottom}{2}{1..160}{bar:hydrophobicity[Red,Gray10]}{}
    %\showfeaturestylename{bbottom}{hydrophob.}

    %\bargraphstretch{3}
    %\featurestylenamescolor{Red}
    %\featurestylenamesrm  \featurestylenamesit

    %\showsequencelogo{top}


    %\showconsensus[ColdHot]{bottom}
    \showconsensus[black]{top}

    %\defconsensus{.}{lower}{upper}
    %\defconsensus{{}}{lower}{upper}
    %\defconsensus{{}}{*}{upper}

    """)
        file_handle.write(features)
        file_handle.write(r"""
    \featurerule{1mm}""")
        if(showlegend):
            file_handle.write(r"""
    \showlegend""")

        align = AlignIO.read(open(aln_fname,'r'), "fasta")
        for a,i in zip(align,range(len(align))):
            # print a.id.replace('|',' | ')
            file_handle.write("""
    \\nameseq{%d}{%s}"""%(i+1,a.id.replace('|',' | ').replace('_','-')))

        file_handle.write(r"""

    \end{texshade}
    """)

def feature_str2dict(featurestring,position='top'):
    """converts string of secondary structure annotation (like in VMD) to our type of dict"""
    #numbering should be 0 based
    #HHHHHHEEEEEBBBBBBCCCCCbTGI
    features=[]
    style=''
    begin=-1
    for i,s in enumerate(featurestring):
        if(s in ['H','G','I']): #helices
            if style!='helix':
                style='helix'
                begin=i
        else:
            if style=='helix':
                style=''
                end=i-1
                features.append({'style':'helix','sel':[begin,end],'position':position})

    for i,s in enumerate(featurestring):
        if(s in ['E','B','b']): #helices
            if style!='-->':
                style='-->'
                begin=i
        else:
            if style=='-->':
                style=''
                end=i-1
                features.append({'style':'-->','sel':[begin,end],'position':position})
    return features




    #prof=cons_prof(alignment)
    #pylab.plot(prof)
if __name__ == '__main__':
    human_h2a_z_core=Seq('SRSQRAGLQFPVGRIHRHLKSRTTSHGRVGATAAVYSAAILEYLTAEVLELAGNASKDLKVKRITPRHLQLAIRGDEELDSLI-KATIAGGGVIPHIHKSLIG')
    xenopus_h2a_core=Seq('TRSSRAGLQFPVGRVHRLLRKGNYAE-RVGAGAPVYLAAVLEYLTAEILELAGNAARDNKKTRIIPRHLQLAVRNDEELNKLLGRVTIAQGGVLPNIQSVLLP')
    
    # human_h2a_z_core=Seq('SRSQRAGLQFPVGRIHRHLKSRTTSHGRVGATAAVYSAAILEYLTAEVLELAGNASKDLKVKRITPRHLQLAIRGDEELDSLIKATIAGGGVIPHIHKSLIG')
    msa=MultipleSeqAlignment([SeqRecord(xenopus_h2a_core,id='H2A',name='H2A'),SeqRecord(human_h2a_z_core,id='H2AZ',name='H2AZ')])
    # features=get_hist_ss_in_aln_for_shade(msa,below=True)
    features=[{'style':'fill:$\\uparrow$','sel':[5,10],'text':'test'}]
    print(features)
    shade_aln2png(msa,filename='default',shading_modes=['charge_functional'], legend=False, features=features,title='',logo=False,hideseqs=False,splitN=20,setends=[],ruler=False,show_seq_names=False,show_seq_length=False)
    



    
            
