import signal
from fastapi import FastAPI, Response
from pydantic import BaseModel
import os
import boto3
import scanpy as sc
import anndata as ad
import pandas as pd
import json

class initializeProjectRequest(BaseModel):
    user:str
    project: str
    datasets: list[str]

class clusteringRequest(BaseModel):
    user: str
    project: str
    resolution: float

class geneRequest(BaseModel):
    user: str
    project: str
    gene_list: list[str]
class annoRequest(BaseModel):
    user: str
    project: str
    annotations: object
# variables in Global context
user_environment = {}
adata = ad.AnnData()
workspace_path = r'/tmp'
s3_plots_dir = ""
s3 = boto3.client("s3")
resolution_global= 0.0

principle_components=50

# set bucket values depending on the environment
def set_user_env():
    global user_environment
    global workspace_path
    # Create a dictionary-like object (similar to R's new.env())
    # Get the environment variable
    environment = os.getenv("ENVIRONMENT", default="dev")
    # Print a message
    print(f"Cellborg Processing Python container running in environment: {environment}")
    # Set bucket names based on the environment
    if environment == "dev":
        DATASET_BUCKET = "cellborgdatasetuploadbucket"
        QC_DATASET_BUCKET = "cellborgqcdatasetbucket"
    else:
        DATASET_BUCKET = f"cellborg-{environment}-datasetupload-bucket"
        QC_DATASET_BUCKET = f"cellborg-{environment}-qcdataset-bucket"

    # Assign bucket names to the user_environment dictionary
    user_environment["dataset_bucket"] = DATASET_BUCKET
    user_environment["qc_dataset_bucket"] = QC_DATASET_BUCKET

    #create temp folder for keeping files in the workspace
    #handle windows
    # if os.name == 'nt':
    #     workspace_path = r'C:\temp' 
    # if not os.path.exists(workspace_path):
    #     os.makedirs(workspace_path)

# load dataset files for QC from S3 bucket to the local file system under /tmp
# TODO: change /tmp to something else as this might cause problems.


def initializeAdata(s3_singlets_path: str, datasets: list[str]):
    global adata
    global user_environment

    samples={datasetId:s3_singlets_path+f"/{datasetId}/singlets.h5ad" for datasetId in datasets} 
    adatas={}
    for sample_id, filepath in samples.items():
        s3.download_file(
            Bucket = user_environment["qc_dataset_bucket"], 
            Key = filepath, 
            Filename= 'singlets.h5ad') 
        sample_adata = sc.read_h5ad("singlets.h5ad")
        adatas[sample_id] = sample_adata
        os.remove('singlets.h5ad')
        print("removed file from local")
        #response = s3.get_object(Bucket= user_environment["qc_dataset_bucket"], Key=filepath)
        #print(f"Pulled {sample_id} from s3: ")
        #file_pulled = response['Body'].read()
        #print(file_pulled)
        #print(type(file_pulled))
        #sample_adata = sc.read_h5ad(file_pulled)
        #sample_adata.var_names_make_unique()
        #adatas[sample_id] = sample_adata
    
    adata = ad.concat(adatas, label="sample")
    adata.obs_names_make_unique()

def normalize():
    global adata
    print("-----normalize begins----")
    # Normalizing to median total counts
    sc.pp.normalize_total(adata)
    # Logarithmize the data
    sc.pp.log1p(adata)
    print ("----normalize completed-----")

def feature_selection():
    global adata
    # ## Feature selection
    # 
    # As a next step, we want to reduce the dimensionality of the dataset and only include the most informative genes. This step is commonly known as feature selection. The scanpy function `pp.highly_variable_genes` annotates highly variable genes by reproducing the implementations of Seurat {cite}`Satija2015`, Cell Ranger {cite}`Zheng2017`, and Seurat v3 {cite}`stuart2019comprehensive` depending on the chosen `flavor`. 
    print("-------- feature selection begins --------")
    sc.pp.highly_variable_genes(adata, n_top_genes=2000)
    sc.pl.highly_variable_genes(adata)

def dimentionality_reduction(s3_path):
    global adata
    global principle_components
    # ## Dimensionality Reduction
    # Reduce the dimensionality of the data by running principal component analysis (PCA), which reveals the main axes of variation and denoises the data.
    print("------- dimentionality reduction begins ------")
    sc.tl.pca(adata)

    #find ideal number of principle components
    df = adata.uns["pca"]['variance_ratio'].cumsum(axis=0)
    if df[-1] < 0.95:
        print("use 50")
    else:
        for num in range(len(df)):
            if df[num] >=0.95:
                principle_components = num+1
                break
                
    # Let us inspect the contribution of single PCs to the total variance in the data. This gives us information about how many PCs we should consider in order to compute the neighborhood relations of cells, e.g. used in the clustering function {func}`~scanpy.tl.leiden` or {func}`~scanpy.tl.tsne`. In our experience, there does not seem to be signifigant downside to overestimating the numer of principal components.
    sc.pl.pca_variance_ratio(adata, n_pcs=principle_components, log=True, save=".png")
    png_file1 = "./figures/pca_variance_ratio.png"
        # Create S3 object key for quality control data
    s3_key = f"{s3_path}/QCpca_variance_ratio.png"
    upload_plot_to_s3(s3_key,png_file1)
    
    # You can also plot the principal components to see if there are any potentially undesired features (e.g. batch, QC metrics) driving signifigant variation in this dataset. In this case, there isn't anything too alarming, but it's a good idea to explore this.
    sc.pl.pca(
        adata,
        color=["pct_counts_mt", "pct_counts_mt"],
        dimensions=[(0, 1), (1, 2)],
        ncols=2,
        size=2,
        save=".png",
    )
    # upload plots to s3
    png_file2 = "./figures/pca.png"
    # Create S3 object key for quality control data
    s3_key = f"{s3_path}/QCpca.png"
    upload_plot_to_s3(s3_key,png_file2)

    print("----- nearest_neighbor_graph begins -----")
    sc.pp.neighbors(adata, n_pcs=principle_components)
    # This graph can then be embedded in two dimensions for visualiztion with UMAP (McInnes et al., 2018):
    sc.tl.umap(adata)
    # We can now visualize the UMAP according to the `sample`. 
    sc.pl.umap(
        adata,
        # Setting a smaller point size to get prevent overlap
        size=2,
        save="1.png"
    )
    targetfile=f"{s3_path}/QCnearest_neighbor.png"
    upload_plot_to_s3(targetfile,"./figures/umap1.png")

def upload_plot_to_s3(s3_key, localfile):
    # Upload the JSON data to S3
    s3.upload_file(localfile, user_environment['qc_dataset_bucket'], s3_key, Callback=print)
    print(f"Uploaded plot png to S3: {s3_key}")

def upload_anndata_to_s3(s3_path, adata):
    # Save the AnnData object to a temporary file
    temp_file = "/tmp/adata.h5ad"
    adata.write(temp_file)
    
    
    bucket_name = user_environment['qc_dataset_bucket']
    s3.upload_file(temp_file, bucket_name, s3_path)

def clustering(s3_path, resolution):
    global adata
    global resolution_global
    global s3_plots_dir

# ## Clustering
# 
# As with Seurat and many other frameworks, we recommend the Leiden graph-clustering method (community detection based on optimizing modularity) {cite}`traag2019louvain`. Note that Leiden clustering directly clusters the neighborhood graph of cells, which we already computed in the previous section.

# Using the igraph implementation and a fixed number of iterations can be significantly faster, especially for larger datasets
    print("------- clustering begins --------")
    sc.tl.leiden(adata, resolution = resolution, flavor="igraph", n_iterations=-1)
    #sc.pl.umap(adata, color=["leiden"], save="2.png")

    # **** creating JSON for clustering ******
    print("test X_umap")
    print(adata.obsm)
    umap_df = pd.DataFrame(
    adata.obsm["X_umap"], columns=["UMAP1", "UMAP2"], index=adata.obs.index
    )
    umap_df["cluster"] = adata.obs["leiden"]


    umap_dict = umap_df.to_dict(orient="index")

    print('creating umap json file...')
    with open("umap_clusters.json", "w") as f:
        json.dump(umap_dict, f)

    #upload umap to s3
    print('uploading umap json to s3...')
    umap_path = f"{s3_path}/UMAP_CLUSTERING&res={int(resolution*100)}.json"
    upload_plot_to_s3(umap_path, 'umap_clusters.json')

    #delete temp file
    print('removing umap json from local...')
    os.remove('umap_clusters.json')

    resolution_global = resolution
    return (
        adata.obs["leiden"].cat.categories
    )

def cell_type_annotation(annotations):
    global adata
    global resolution_global
    
    adata.obs["cell_type_lvl1"] = adata.obs["leiden"].map(annotations)
#----- main -------
app = FastAPI()
#------------------- Processing & Annotations-----------------

@app.post("/init_endpoint", status_code=200)
async def initialize_project(initReq: initializeProjectRequest):
    try:
        global adata
        global user_environment
        global principle_components

        s3_path = f"{initReq.user}/{initReq.project}"
        set_user_env()
        initializeAdata(s3_path, initReq.datasets)
        print('Successfully concatonated all datasets')
        normalize()
        print('Successfully normalized concatonated dataset')
        feature_selection()
        print('Successfully selected features')
        dimentionality_reduction(s3_path)
        print('Successfully reduced dimensions')
        
        #create project_values.json in proj dir
        numpcs = {
            "num_pcs":principle_components,
            "gene_list": adata.var_names.to_list()
        }
        with open("project_values.json", "w") as outputfile:
            json.dump(numpcs, outputfile)
        print("created project_values.json file")

        input_var_path = f"{initReq.user}/{initReq.project}/project_values.json"
        upload_plot_to_s3(input_var_path, 'project_values.json')
        print("uploaded project_values.json file to s3")

        os.remove("project_values.json")
        print("project_values.json file successfully deleted")

        return{
            "success": True,
            "message": "ProcAnno successfully initialized"
        }

    except Exception as err:
        print('ERROR: ',str(err))
        return {
            "success": False,
            "message": str(err)
        }

@app.post("/clustering", status_code=200)
async def do_clustering(clustReq: clusteringRequest):
    try:
        s3_path = f"{clustReq.user}/{clustReq.project}"
        clusters = clustering(s3_path, clustReq.resolution)

        #download old project_values file from s3
        s3.download_file(
            Bucket = user_environment["qc_dataset_bucket"], 
            Key = s3_path+'/project_values.json', 
            Filename= 'project_values.json') 

        #add clustering resolution
        with open('project_values.json', "r+") as f:
            data = json.load(f)
            data['clust_resolution'] = clustReq.resolution
            #write over file
            f.seek(0)
            json.dump(data, f)
            f.truncate()
        print("Added clustering resolution to project values")
        
        #upload updated file to s3
        upload_plot_to_s3(f"{s3_path}/project_values.json",'project_values.json')
        print("uploaded new project values to s3")
        
        #remove temp file locally
        os.remove("project_values.json")

        clusters = clusters.to_list()
        print(f"cluster: {type(clusters)}", clusters)
        return {
            "success": True,
            "message": "Clustering successfully finished",
            "clusters": clusters
        }
    except Exception as err:
        print('ERROR: ',str(err))
        return {
            "success": False,
            "message": str(err)
        }

@app.post("/gene_expression", status_code = 200)
async def gene_expression(geneReq: geneRequest):
    global adata
    try:
        gene_list = geneReq.gene_list
        print('----Starting Gene Expression----')
        gene_mask = adata.var_names.isin(gene_list)
        selected_genes = adata.var_names[gene_mask]
        expression_df = pd.DataFrame(
        adata[:, gene_mask].X.toarray(),  # Convert to dense matrix if needed
        columns=selected_genes,
        index=adata.obs.index
        )
        expression_df["cluster"] = adata.obs["leiden"]
        umap_df = pd.DataFrame(
        adata.obsm["X_umap"], columns=["UMAP1", "UMAP2"], index=adata.obs.index
        )
    
        combined_df = pd.concat([umap_df, expression_df], axis=1)
        combined_dict = combined_df.to_dict(orient="index")

        print("saving gene expression to local...")
        with open("gene_expression_per_cell_with_clusters.json", "w") as f:
            json.dump(combined_dict, f)
        
        s3_path = f"{geneReq.user}/{geneReq.project}"

        upload_plot_to_s3(f"{s3_path}/gene_expression.json","gene_expression_per_cell_with_clusters.json")
        print("uploaded gene expression json to s3")

        print("----removing temp file")
        os.remove("gene_expression_per_cell_with_clusters.json")

        return{
            "success": True,
            "message": "Successfully conducted gene expression"
        }
    except Exception as err:
        print('ERROR: ',str(err))
        return {
            "success": False,
            "message": str(err)
        }

@app.post("/annotations", status_code = 200)
async def annotations(annotateRequest: annoRequest):
    global adata
    #try:
    print("------Starting annotations------")
    #annotations_dict = {str(i): annotateRequest.annotations[i] for i in range(len(annotateRequest.annotations))}
    cell_type_annotation(annotateRequest.annotations)
    #used to verify that annotations did work
    print("creating test png")
    sc.pl.umap(
    adata,
    color=["cell_type_lvl1"],
    legend_loc="on data",
    save = "annotations_test.png"
    )
    print("uploading test png")
    upload_plot_to_s3(f"{annotateRequest.user}/{annotateRequest.project}/annotations_test.png", "./figures/umapannotations_test.png")
    return{
        "success":True,
        "message":"Annotatons successfully completed"
    }
    #except Exception as err:
    #    print('ERROR: ',str(err))
    #    return{
    #        "success":False,
    #        "message": str(err)
    #    }

@app.post("/shutdown")
async def shutdown(user: str, project: str):
    global adata
    try:
        if user and project:
            # Assuming adata is your AnnData object and you have a function to upload it to S3
            upload_anndata_to_s3(f"{user}/{project}/adata.h5ad", adata)
            print("AnnData object uploaded to S3 successfully.")
        else:
            print('User and project not provided. AnnData object not uploaded to S3.')
    except Exception as e:
        print(f"Failed to upload AnnData object to S3: {str(e)}")
    
    os.kill(os.getpid(), signal.SIGTERM)
    return Response(status_code=200, content='Server shutting down...')

@app.get("/health", status_code = 200)
async def health():
    
    return {"status": "ok"}