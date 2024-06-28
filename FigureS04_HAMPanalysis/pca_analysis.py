import pandas as pd
import numpy as np
import os
from prody import parsePDB, AtomGroup
from Bio import SeqIO
import warnings
warnings.filterwarnings("ignore")
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt


def get_msa_mask(msa_file):
    # Read the MSA output file from foldseek
    def masking_array(binary_array):
        # Create a mask that identifies columns containing at least one 0
        mask = np.any(binary_array == 0, axis=0)
        # Set all elements in columns containing at least one 0 to 0
        binary_array[:, mask] = 0
        return binary_array
        
    fasta_sequences = list(SeqIO.parse(open(msa_file),'fasta'))
    for i, record in enumerate(fasta_sequences):
        seq = np.array(list(record.seq))
        id_ = record.id
        
        # Initialize the masks array that indicates where the "-" is in the alignment
        if i==0:
            masks = np.zeros((len(fasta_sequences), len(seq)))
        if len(seq) > len(masks[i]):
            n = int(len(seq)-len(masks[i]))
            masks = np.pad(masks, ((0, 0), (0, n)), mode='constant', constant_values=0)
        if len(seq) < len(masks[i]):
            masks[i][0:len(seq)] = np.array([0 if "-"==element else 1 for element in seq])
            continue        
        masks[i] = np.array([0 if "-"==element else 1 for element in seq])
    return masking_array(masks)

def get_rototranslation(line):
    T = np.array([line[2:5]])
    R = np.array([line[5:8], line[8:11],line[11:14]])
    M = np.eye(4) 
    M[:3, :3] = R
    M[:3, 3] = T
    return M

def read_aln(pdbs_folder, filename):
    coords = {}
    df = pd.read_csv(filename, header=None, sep=" ")
    df[["Query","Target"]]= df[0].str.split("\t", expand=True)
    queries = df.Query.unique()
    for query in queries:
        tmp = df[df.Query==query]
        new = {}
        for _,line in tmp.iterrows():
            # Get the transformation matrix of the target protein
            transformation_matrix = get_rototranslation(line)
            # Identify the pdb file of the protein
            pdb_file = os.path.join(pdbs_folder,line.Target+".pdb")
            # Transform the target protein and return it as a dictionary element
            new[line.Target] = run_rototranslations(pdb_file, "A", transformation_matrix)
        coords[query.rsplit("_", 1)[1]] = new
    # Return the final dictionary
    return coords

def run_rototranslations(pdb_file, chain, transformation_matrix):
    
    # Read protein file
    protein = parsePDB(pdb_file, chain=chain, subset="ca")

    # Get the coordinates of the selected atoms
    coords = protein.getCoords()

    # Convert the coordinates to homogeneous coordinates by adding a column of ones
    homogeneous_coords = np.hstack((coords, np.ones((coords.shape[0], 1))))

    # Perform the transformation
    transformed_coords = np.dot(homogeneous_coords, transformation_matrix)

    # Convert back to Cartesian coordinates by removing the homogeneous coordinate
    transformed_coords_cartesian = transformed_coords[:, :3]
    
    # Return the transformed coordinates in a pandas dataframe
    return pd.DataFrame(transformed_coords_cartesian)
    

def clean_up_coordset(mask, coord_set):
    min_ = np.where(mask==1)[0][0]
    max_ = np.where(mask==1)[0][-1]
    flags = np.zeros(len(coord_set))
    counter = 0
    for i in range(min_,max_+1):
        if mask[i]==1:
            flags[counter]=1.0
            counter += 1
    coord_set["Flag"] = flags
    return coord_set[coord_set.Flag==1]

def prepare_coordset_for_pca(masks, pdbs):
    for i,key in enumerate(pdbs.keys()):
          coords = clean_up_coordset(masks[i+1], pdbs[key])
          coords.drop(["Flag"], axis=1, inplace=True)
          if i==0: pca = coords.to_numpy().flatten()
          else: pca = np.vstack([pca, coords.to_numpy().flatten()])
    return pca

fig, axs = plt.subplots(1, 1, figsize=(7, 5))  # Adjust figsize as needed
plt.rcParams.update({'font.family': 'Arial'})
cm = plt.cm.get_cmap('viridis')
pdbs = read_aln("struct_in/multimers", "Foldseek/Multimer/Alignment/aln_tmscore.tsv")

def map_ids(data):
    # Get unique elements and their indices
    unique_elements, unique_indices = np.unique(data, return_inverse=True)
    # Generate unique integer values
    unique_integers = np.arange(len(unique_elements))
    # Map original data to unique integers
    return unique_integers[unique_indices]


for i,(axs,chain) in enumerate(zip([axs], pdbs.keys())):
    chain_set = prepare_coordset_for_pca(get_msa_mask(f"Foldseek/Multimer/MSA/msa_output/{i}.a3m"), pdbs[chain])
    # Perform PCA with the desired number of components
    num_components = 2  # You can change this to the number of components you want
    pca = PCA(n_components=num_components)
    pca_result = pca.fit_transform(chain_set)

    axs.scatter(pca_result[0, 0], pca_result[0, 1], c="red", marker="x", s=150, linewidths=3, label="AlphaFold HAMP")
    mapped_ids = map_ids([i.rsplit("_", 1)[0] for i in pdbs[chain].keys()])
    map_ = axs.scatter(pca_result[1:, 0], pca_result[1:, 1], c=mapped_ids[1:], cmap=cm, s=150)
    axs.set_xlabel("PC1 [{:.2f}]".format(pca.explained_variance_ratio_[0]), fontsize=17)
    axs.set_ylabel("PC2 [{:.2f}]".format(pca.explained_variance_ratio_[1]), fontsize=17)
    axs.spines[['right', 'top']].set_visible(False)
    axs.spines["bottom"].set_linewidth(1)
    axs.spines["left"].set_linewidth(1)
    axs.tick_params(length=5, width=1, labelsize=15)

plt.legend(loc="upper right", ncol=1, fontsize=15)

cbar = plt.colorbar(map_, extend='both')
cbar.ax.tick_params(labelsize=12)
bar_labels = []
for i in [x for x in pdbs[chain].keys()][1:]:
    pdb_id = i.rsplit("_", 1)[0]
    if pdb_id not in bar_labels:
        bar_labels.append(pdb_id)
cbar.ax.locator_params(nbins=len(bar_labels))
cbar.ax.set_yticklabels(bar_labels)

fig.tight_layout()
plt.savefig("FigureS04_HAMPanalysis.png")
