import pandas as pd
import openai
import sys


def get_cluster_cell_types(input_filepath, input_metadata):
    # Read the file and convert the dataframe to a string
    df = pd.read_csv(input_filepath, sep="\t")
    df = df[["gene", "p_val_adj", "avg_log2FC", "cluster"]]

    # Convert the dataframe to a string
    df_str = df.to_string()

    try:
        openai.api_key = "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"
    except KeyError:
        sys.stderr.write(
            """
            Error: Missing OPENAI_API_KEY environment variable
            """
        )
        exit(1)

    response = openai.ChatCompletion.create(
        model="gpt-4",
        messages=[
            {
                "role": "system",
                "content": "You are a bioinformatics expert. Your role is to annotate clusters by cell type. The input data is a tab-separated table with information about the most highly expressed genes in a given cluster. You are to return predicted cell types of a given cluster in the format of a character vector. Clusters can be made up of multiple cell types or unknown cell types. Use your knowledge about marker genes in specific tissue. Only return cell types specific to the tissue of interest.",
            },
            {
                "role": "user",
                "content": f"The gene expression data belongs to {input_metadata}.",
            },
            {
                "role": "assistant",
                "content": "Ok I am ready for the input table. I will return the predicted cell type of each cluster based upon the marker gene, p_val_adj (p-value adjusted) column and the avg_log2fc (average log2 fold change) column. I will only return the cluster number and the predicted cell type.",
            },
            {
                "role": "user",
                "content": f"Recall which tissue the data comes from: {input_metadata}. Be precise with your annotation. Only return cluster and cell type prediction. Here is the data: \n{df_str}",
            },
        ],
    )

    return response
