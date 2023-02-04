import shutil
from io import StringIO
import streamlit as st
import plotly.express as px
import streamlit as st
from Bio.PDB import PDBParser, PDBList, MMCIFParser
import numpy as np
import math
import pandas as pd
import os
import streamlit.components.v1 as components
########################################################################################################################
# Making WebPage Wide for App and Setting Name of app to be displayed on tab in browser

st.set_page_config(
    page_title="Ramachandran",
    page_icon="favicon2.png",
    layout="wide",
    #initial_sidebar_state="expanded",
    #menu_items={'About': "# This App has been developed by Junaid"}
)

##############################    TITLE   #####################################################################

st.markdown("<center><br><h1 style='background-color:#330000; color:#ffdb99'>Torsion Angles and Ramachandran Plot<br></h1></center>",unsafe_allow_html=True)


###########################    DELETING FILES FROM LOCAL DIRECTORY UPON REFRESH    #############################

folder = 'pdb_files/'
try:
    for file_name in os.listdir(folder):
        file_path = os.path.join(folder, file_name)
        try:
            if os.path.isfile(file_path) or os.path.islink(file_path):
                os.unlink(file_path)
            elif os.path.isdir(file_path):
                shutil.rmtree(file_path)
        except Exception as e:
            print('Failed to delete %s. Reason: %s' % (file_path, e))
except:
    pass

#######################     TORSION ANGLE CALCULATION      #############################################################
@st.cache
def torsion_angles(p1, p2, p3, p4):
    """Calculate the torsion angle between plane p1-p2-p3 and plane p2-p3-p4.
    p1,p2,p3,p4 are three dimensional coordinates"""
    b1 = p2 - p1
    b2 = p3 - p2
    b3 = p4 - p3
    if (np.cross(b1, b2) == (0.0, 0.0, 0.0)).all() or (np.cross(b2, b3) == (0.0, 0.0, 0.0)).all():
        return 180.0
    else:
        n1 = np.cross(b1, b2) / np.linalg.norm(np.cross(b1, b2))
        n2 = np.cross(b2, b3) / np.linalg.norm(np.cross(b2, b3))

        m1 = np.cross(n1, b2 / np.linalg.norm(b2))
        x = np.dot(n1, n2)
        y = np.dot(m1, n2)
        return -math.atan2(y, x) / math.pi * 180.0

#####################       MAKE DATAFRAME OF TORSION ANGLES    ########################################################
@st.cache
def make_dataframe(fname):
    try:
        file = files_dict.get(fname)
        file = StringIO(file.getvalue().decode("utf-8"))
    except:
        path = 'pdb_files/'
        file = path + f'{fname}'
    finally:
        pass

    try:
        parser = MMCIFParser(QUIET=True)
        structure = parser.get_structure(f'{fname}', file)
    except:
        parser = PDBParser(PERMISSIVE=True,QUIET=True)
        structure = parser.get_structure(f'{fname}', file)
        
    list_n_coord = []
    list_ca_coord = []
    list_c_coord = []
    list_res = []
    list_res_chain = []
    list_res_seq = []
    list_discontinuous_seq_f = []
    list_discontinuous_seq_r = []
    list_noinfo_seq_f = []
    list_noinfo_seq_r = []
    list_aa = ['GLY', 'ALA', 'VAL', 'LEU', 'ILE', 'PRO', 'PHE', 'TYR', 'TRP', 'SER',
               'THR', 'CYS', 'MET', 'ASN', 'GLN', 'ASP', 'GLU', 'LYS', 'ARG',
               'HIS']  # Build a table of 20 common amino acids for screening
    models = structure.get_list()
    #"""
    #By traversing the various subclasses of the parser，
    #1) Record the three-letter name of each AA;
    #2) Record all N, CA, C space coordinate arrays required to calculate the dihedral angle and store as a list.
    #"""
    not_continue_idx = []
    no_info_idx = []
    idx = 0
    real_idx = 0
    need_atom = ['N', 'CA', 'C']
    for model in models:
        chains = model.get_list()
        for chain in chains:
            residues = chain.get_list()
            for residue in residues:
                res_name = residue.get_resname()
                if res_name in list_aa:  # Change the previous "res_name != 'UNK'" to determine whether it is in common amino acids
                    list_res.append(res_name)
                    atoms = residue.get_list()
                    list_res_chain.append(
                        residue.get_full_id()[2])  # record the chain in which the amino acid is located
                    list_res_seq.append(residue.get_id()[1])  # Record the order of appearance of amino acids

                    all_atom = []
                    for atom in atoms:
                        atom_name = atom.get_name()
                        all_atom.append(atom_name)
                    if all([i in all_atom for i in need_atom]):
                        for atom in atoms:
                            atom_name = atom.get_name()
                            if atom_name == 'N':
                                atom_coord = atom.get_coord()
                                list_n_coord.append(atom_coord)
                            elif atom_name == 'CA':
                                atom_coord = atom.get_coord()
                                list_ca_coord.append(atom_coord)
                            elif atom_name == 'C':
                                atom_coord = atom.get_coord()
                                list_c_coord.append(atom_coord)
                        if len(list_res_seq) > 1:
                            if residue.get_id()[1] != (
                                    list_res_seq[-2] + 1):  # Determine if amino acids are consecutive
                                # Record whether it is disconnected from the previous amino acid
                                not_continue_idx.append(idx)
                                list_discontinuous_seq_f.append(list_res_seq[-2])  # Record the sequence number of the previous amino acid at the discontinuity
                                list_discontinuous_seq_r.append(residue.get_id()[1])  # Record the sequence number of the next amino acid at the discontinuity
                        idx += 1
                    else:
                        no_info_idx.append(real_idx)
                        # list_noinfo_seq_f.append(list_res_seq[-2]) # 记录不连续处前一氨基酸序号
                        # list_noinfo_seq_r.append(residue.get_id()[1]) # 记录不连续处后一氨基酸序
                    real_idx += 1

                    # calculate two dihedral angles
    n_numb = len(list_n_coord)
    ca_numb = len(list_ca_coord)
    c_numb = len(list_c_coord)

    list_psi = []
    list_phi = ['N/A']  # phi angle is not applicable for the first amino acid
    list_omega = ['N/A']  # omega angle does not apply for the first amino acid
    for i in range(0, (len(list_n_coord) - 1)):
        n1 = list_n_coord[i]
        ca1 = list_ca_coord[i]
        c1 = list_c_coord[i]
        n2 = list_n_coord[i + 1]
        ca2 = list_ca_coord[i + 1]
        c2 = list_c_coord[i + 1]
        psi = torsion_angles(n1, ca1, c1, n2)
        phi = torsion_angles(c1, n2, ca2, c2)
        omega = torsion_angles(ca1, c2, n2, ca2)
        list_psi.append(psi)
        list_phi.append(phi)
        list_omega.append(omega)
    list_psi.append('N/A')  # psi angle is not applicable for the first amino acid

    for i in no_info_idx:
        list_res[i] = 'noinf'
        list_res_chain[i] = 'noinf'
        list_res_seq[i] = 'noinf'
    for x in list_res_chain:
        while list_res_chain.count('noinf') > 0:
            list_res_chain.remove('noinf')
    for x in list_res_seq:
        while list_res_seq.count('noinf') > 0:
            list_res_seq.remove('noinf')
    for x in list_res:
        while list_res.count('noinf') > 0:
            list_res.remove('noinf')

    # Invalidate the phi angle of the last amino acid at all breaks, and invalidate the phi angle and omega angle of the first amino acid


    for i in not_continue_idx:
        list_phi[i] = 'N/A'
        list_psi[i - 1] = 'N/A'
        list_omega[i] = 'N/A'

    Dict = {'chain': list_res_chain, 'resSeq': list_res_seq, 'Amino Acids': list_res,
            'psi': list_psi, 'phi': list_phi, 'omega': list_omega}
    Dict['file'] = fname

    new_table = pd.DataFrame(Dict, columns=['chain', 'resSeq', 'Amino Acids', 'psi', 'phi', 'omega', 'file'])

    return new_table

############################        MAKE Plot           ################################################################

def make_plot(PDB_file,Torsion_x,Torsion_y):

    df = table[table['file'] == PDB_file]
    df = df[df.psi != 'N/A']
    df = df[df.phi != 'N/A']
    df = df[df.omega != 'N/A']
    df = df[df['Amino Acids'] != 'N/A']

    Gly = df[df['Amino Acids'] == 'GLY']
    Pro = df[df['Amino Acids'] == 'PRO']
    pro_index = df[df['Amino Acids'] == 'PRO'].index.tolist()
    
    try:
        pre_pro_index = [i - 1 for i in pro_index]
        Pre_pro = df.iloc[pre_pro_index]
    except:
        pass

    if Aa == 'All':
        df = df
    elif Aa == 'Glysine':
        df = Gly
    elif Aa == 'Proline':
        df = Pro
    else:
        df = Pre_pro


    df['psi'] = df['psi'].astype(float)
    df['phi'] = df['phi'].astype(float)
    df['omega'] = df['omega'].astype(float)
    df['Amino Acids'] = df['Amino Acids'].astype(str)

    ############################        PLOTlY       ###############################

    fig2 = px.density_contour(df, x=Torsion_x, y=Torsion_y, range_x=[-180,180],range_y=[-180,180],title=Aa + ' AA')  #color_discrete_sequence=px.colors.qualitative.Light24
    fig2.update_traces(contours_coloring="fill", contours_showlabels = False, colorscale='blackbody')   #colorscale='blackbody'
    fig2.update_traces(showscale=False)

    with maker_col:
        marker_size = st.slider(label='MarkerSize', min_value=0.5, max_value=9.0, value=5.0, step=0.2)
    with opacity_col:
        opacity = st.slider(label='Opacity', min_value=0.0, max_value=1.0, value=0.9, step=0.1)


    fig2.add_trace(px.scatter(df, x=Torsion_x, y=Torsion_y,hover_data=['Amino Acids'],opacity=opacity).data[0])
    fig2.update_traces(marker=dict(size=marker_size,
                                   color='#00ff00',
                                   line=dict(width=1, color='#003311')),
                      selector=dict(mode='markers'))
    fig2.add_hline(y=0,line_width=0.6, line_color="white")
    fig2.add_vline(x=0, line_width=0.6, line_color="white")
    fig2.update_layout(xaxis_range=[-180,180],
                       yaxis_range=[-180,180],
                       paper_bgcolor="rgba(0,0,0,0)",
                       plot_bgcolor="black",
                       title_font={'family': "Ariel",
                                   'color': '#ffbd45',
                                   'size': 39},
                       title_xanchor='auto',
                       title_yanchor='top',
                       margin_autoexpand=True,autosize=True,
                       height= 650, #width = 700
                       )
    fig2.update_xaxes(autorange=True,dtick=45,mirror=True,showline=True, automargin=True)
    fig2.update_yaxes(autorange=False,showgrid=False,dtick=45,mirror=True,showline=True,automargin=True)


    result_col2.title(':orange[Torsion Angles]')

    #return df, fig2
    #return st.dataframe(df),st.plotly_chart(fig2,theme=None, config = {'displayModeBar': False})
    return result_col2.dataframe(df), result_col1.plotly_chart(fig2, theme=None, config={'displayModeBar': False,'scrollZoom': True})

##########################      Downoad Input PDB IDs Files    #########################################################
#@st.cache
def download_PDB(pdb_id_list):
    pdb_id_list = pdb_id_list.split(',')
    pdb = PDBList()
    for pdb_id in pdb_id_list:
        try:
            pdb.retrieve_pdb_file(pdb_id,pdir='pdb_files/')
        except:
            try:
                pdb.retrieve_pdb_file(pdb_id, pdir='pdb_files/',file_format='pdb')
            except:
                st.warning('Something is wrong either with format or ID so please re-validate PDB ID or retry')
        else:
            pass
###############################     Creating File Records       ######################################################

def file_records(fnames,files_dict,uploaded_files):
    #fnames = []
    #files_dict = {}
    if len(uploaded_files) != 0:
        fnames.clear()
        try:
            #fnames = []
            #files_dict = {}
            for uploaded_file in uploaded_files:
                fnames.append(uploaded_file.name)
            for i in range(len(fnames)):
                files_dict[fnames[i]] = uploaded_files[i]
        except:
            pass
    else:
        # fnames = []
        fnames.clear()
        path = 'pdb_files/'
        files = os.listdir(path)
        #st.warning(files)

        for file in files:
            if (file[-4:] == '.pdb') or (file[-4:] == '.ent') or (file[-4:] == '.cif'):
                fnames.append(file)
                # fnames.append(file[0: -4])
    return fnames,files_dict

##############################     UPLOAD FILE WIGET     #################################

uploaded_files = st.file_uploader(label='Upload PDB File', type=['.pdb','.ent','cif'], accept_multiple_files=True)


########################    CONTAINER 1 AND THE COLUMNS      ############################

with st.container():
    col1,col2 = st.columns([2,1],gap='small')
    with col1:
        pdb_id_list = st.text_input("PDB IDs", placeholder="Enter PDB IDs example: 6VMK, 4WZJ,1A0K,3CMA,6JMD, 7SC2, 6LP6")
        download_PDB(pdb_id_list=pdb_id_list)       # Download PDB Files
    with col2:
        fnames = []
        files_dict = {}
        fnames, files_dict = file_records(fnames=fnames, files_dict=files_dict, uploaded_files=uploaded_files)
        fnames.sort()
        PDB_file = st.selectbox('PDB files', options=fnames, index=0)

####################      CREATING DATAFRAME OF TORSION ANGLES      ####################

table = pd.DataFrame(columns=['chain', 'resSeq', 'Amino Acids', 'psi', 'phi', 'omega', 'file'])
fnames.sort()
value_error = []
index_error = []
if len(fnames) != 0:
    for fname in fnames:
        try:
            new_table = make_dataframe(fname)
            #table = table.append(new_table)
            table = pd.concat([table,new_table])
        except:
            st.warning('Something is wrong. Issue may be with format in PDB file.')
            pass
else:
    pass
    #st.warning('No Uploaded File: Either Upload PDB file or Enter PDB ID')

##########################    CONTAINER 3 FOR COL1,2,3    ###################################

with st.container():
    fnames.sort()
    col1,col2,col3 =  st.columns(3,gap='small')
    with col1:
        Aa = st.selectbox('Amino Acids', options=['All', 'Glysine', 'Proline', 'Pre-Proline'], index=0)
    with col2:
        Torsion_x = st.selectbox('Torsions x', options=['phi', 'psi', 'omega'], index=0)
    with col3:
        Torsion_y = st.selectbox('Torsions y', options=['psi', 'phi'], index=0)


##################################### DISPLAYING RAMCHANDAN PLOT AND DATAFRAME #################################

###### DFINING COLUMNS FOR MARKER AND RESULTS AND WILL BE USED INSIDE MAKE PLOT FUNCTION ###################

maker_col, opacity_col = st.columns(2)
result_col1, result_col2 = st.columns([1.5,1])

#col1, col2 = st.columns(2)
#with col2: st.dataframe(df)

if len(fnames) != 0:

    df, figure = make_plot(PDB_file, Torsion_x, Torsion_y)
    #with col1:
    #    st.plotly_chart(figure, theme=None, config={'displayModeBar': False})
    #with col2:
    #    st.title(':orange[Torsion Angles Only For Select File]')
    #    st.dataframe(df)

    ##############      DOWNLOAD DATAFRAME AS CSV    ##################

    @st.cache
    def convert_df(table):
       return table.to_csv(index=False).encode('utf-8')
    data_csv = convert_df(table)
    result_col2.download_button(label= "Download Torsions",data = data_csv,file_name= "torsion_for_all_files.csv",
                   mime= "text/csv", key='download-csv')
else:
    st.warning('Either Upload Protein PDB File/s or Enter PDB IDs ')


##############   Hide hamburger menu    ##############

hide_menu_style = """
        <style>
        #MainMenu {visibility: hidden;}
        </style>
        """
st.markdown(hide_menu_style, unsafe_allow_html=True)

#################### Hide default footer #############

hide_default_footer = """
        <style>
        footer {visibility: hidden;}
        </style>
       """
st.markdown(hide_default_footer, unsafe_allow_html=True)

################################    THE END   THE END         ###########################################################
